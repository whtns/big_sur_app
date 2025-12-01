import os
import shutil
import time
import pickle
from datetime import datetime

from filelock import Timeout, FileLock
import zarr
from numcodecs import Blosc

import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc
import anndata as ad
# Prefer public anndata I/O. Older code used anndata._io.zarr internals
# (read_dataframe/read_attribute/write_attribute) which are not stable.
# Use the public API (`ad.read_zarr`, `AnnData.write_zarr`) instead.
_ANNDATA_ZARR_INTERNALS = False

from status.status_functions import *
from plotting.multi_color_scale import MultiColorScale

from tasks.tasks import write_dense

def safe_write_zarr(adata_obj, zarr_cache_dir):
    """Safely write AnnData to a zarr store.

    Strategy:
    - Acquire a FileLock on the target store path used by both foreground
      code and Celery tasks.
    - Write the AnnData to a temporary zarr directory (atomic replacement).
    - Rename the temp directory over the target to avoid partial-group races.
    - If that fails, attempt best-effort cleanup of known conflicting groups
      (e.g. `X_dense`, `layers_dense`) and retry once.
    """
    lock_filename = str(zarr_cache_dir) + ".lock"
    lock = FileLock(lock_filename, timeout=lock_timeout)

    tmp_dir = None
    try:
        with lock:
            # write to a temp zarr directory then atomically swap
            ts = int(time.time() * 1000)
            pid = os.getpid()
            tmp_dir = f"{zarr_cache_dir}.tmp.{pid}.{ts}"
            try:
                # ensure any previous tmp is removed
                if os.path.exists(tmp_dir):
                    shutil.rmtree(tmp_dir)
            except Exception:
                pass

            try:
                adata_obj.write_zarr(tmp_dir)
            except Exception as e:
                print(f"[WARN] adata.write_zarr to tmp failed: {e}")
                # fall through to group cleanup fallback below
                raise

            # At this point tmp_dir is fully written. Replace the target.
            try:
                if os.path.exists(zarr_cache_dir):
                    # remove old store atomically
                    backup = f"{zarr_cache_dir}.bak.{pid}.{ts}"
                    try:
                        if os.path.exists(backup):
                            shutil.rmtree(backup)
                    except Exception:
                        pass
                    try:
                        os.rename(zarr_cache_dir, backup)
                    except Exception:
                        # fallback to remove
                        shutil.rmtree(zarr_cache_dir)
                os.rename(tmp_dir, zarr_cache_dir)
                # remove backup if present
                try:
                    if 'backup' in locals() and os.path.exists(backup):
                        shutil.rmtree(backup)
                except Exception:
                    pass
                return
            except Exception as e:
                print(f"[WARN] atomic replace of zarr store failed: {e}")
                # try best-effort cleanup of conflicting groups below

    except Exception:
        # lock context or tmp write failed — fall through to group cleanup
        pass

    # Fallback: try removing conflicting helper groups then retry write
    try:
        store_store = zarr.DirectoryStore(zarr_cache_dir)
        grp = zarr.open_group(store=store_store, mode='a')
        removed = False
        for gname in ["X_dense", "layers_dense"]:
            if gname in grp:
                try:
                    del grp[gname]
                    print(f"[DEBUG] removed conflicting group {gname} from zarr store")
                    removed = True
                except Exception as e2:
                    print(f"[DEBUG] failed removing {gname}: {e2}")
        if removed:
            try:
                # try writing directly to target now
                adata_obj.write_zarr(zarr_cache_dir)
                return
            except Exception as e3:
                print(f"[ERROR] adata.write_zarr retry after cleanup failed: {e3}")
    except Exception as e4:
        print(f"[DEBUG] opening zarr store for cleanup failed: {e4}")

    # If we reach here, surface an exception
    raise RuntimeError(f"safe_write_zarr: failed to write zarr to {zarr_cache_dir}")

save_analysis_path = "/Library/WebServer/Documents/BigSuR/cache/"
selected_datasets_path = "/Library/WebServer/Documents/BigSuR/selected_datasets/"
user_dataset_path  = "/Library/WebServer/Documents/BigSuR/user_datasets/"

lock_timeout = 60

use_zarr = True

### the actual helper functions
### TODO: break all of this up into specific modules

def generate_adata_from_10X(session_ID, data_type="10X_mtx"):
    data_dir = save_analysis_path + str(session_ID) + "/raw_data/"
    if (data_type == "10X_mtx"):
        adata = sc.read_10x_mtx(data_dir, cache=False)
    elif (data_type == "10X_h5"):
        adata = sc.read_10x_h5(data_dir + "data.h5ad")
    else:
        print("[ERROR] data type not recognized - returning None")
        return None

    cache_adata(session_ID, adata)
    return adata

def load_selected_dataset(session_ID, dataset_key):
    dataset_dict = {
    "00000": "pbmc3k_raw",
    "00001": "pbmc3k_processed",
    }

    filename = dataset_dict[dataset_key]
    if (filename is None):
        return None
    else:
        filename = selected_datasets_path + filename
    print(filename)
    adata = sc.read_h5ad(filename + ".h5ad")

    state = {"filename": str(dataset_dict[dataset_key]),
             "# cells/obs": len(adata.obs.index),
             "# genes/var": len(adata.var.index),
             "# counts": int(np.sum(adata.X))}
    cache_state(session_ID, state)

    adata = cache_adata(session_ID, adata)
    return adata

import ipdb
def cache_adata(session_ID, adata=None, group=None,
                store_dir=None, store_name=None):
    if ((store_dir is None) or (store_name is None)):
        save_dir = save_analysis_path  + str(session_ID) + "/"
        filename = save_dir + "adata_cache"
        chunk_factors = [150, 3] #faster, hot storage
    else:
        save_dir = store_dir
        filename = save_dir + store_name
        chunk_factors = [3, 3] #slower, cold storage

    
    if not (os.path.isdir(save_dir)):
        try:
            print("[DEBUG] making directory:" + str(save_dir))
            os.mkdir(save_dir)
        except:
            return None
    
    lock_filename = (save_analysis_path  + str(session_ID) 
                     + "/" + "adata.lock")
    lock = FileLock(lock_filename, timeout=lock_timeout)

    compressor = Blosc(cname='blosclz', clevel=3, 
                       shuffle=Blosc.SHUFFLE)
    zarr_cache_dir = filename  + ".zarr"
    attribute_groups = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "uns", "raw"]
    extra_attribute_groups = ["X_dense", "layers_dense"]

    if (adata is None): # then -> read it from the store
        if (os.path.exists(zarr_cache_dir) is True):
            # Use public anndata read API to load the dataset. ad.read_zarr
            # returns a full AnnData; then return either the requested
            # attribute group or the whole AnnData.
            try:
                adata = ad.read_zarr(zarr_cache_dir)
            except Exception as e:
                # Provide diagnostic information to help debug the zarr store
                print(f"[ERROR] ad.read_zarr failed: {e}")
                try:
                    store_store = zarr.DirectoryStore(zarr_cache_dir)
                    store = zarr.open_group(store=store_store, mode='r')
                    try:
                        keys = list(store.group_keys())
                    except Exception:
                        try:
                            keys = list(store.keys())
                        except Exception:
                            keys = None
                    print(f"[DEBUG] zarr store keys: {keys}")
                except Exception as e2:
                    print(f"[DEBUG] opening zarr store failed: {e2}")
                return None

            if (group in attribute_groups):
                mapping = {
                    'obs': adata.obs,
                    'var': adata.var,
                    'obsm': adata.obsm,
                    'varm': adata.varm,
                    'obsp': adata.obsp,
                    'varp': adata.varp,
                    'layers': adata.layers,
                    'X': adata.X,
                    'uns': adata.uns,
                    'raw': adata.raw
                }
                return mapping.get(group, None)
            else:
                return adata
    else: # then -> update the state dictionary and write adata to the store
        if (group is None):
            cache_state(session_ID, key="# cells/obs", val=len(adata.obs.index))
            cache_state(session_ID, key="# genes/var", val=len(adata.var.index))
            if ("total_counts" in adata.obs):
                cache_state(session_ID, key="# counts", val=int(np.sum(adata.obs["total_counts"])))
            else:
                cache_state(session_ID, key="# counts", val=int(np.sum(adata.X)))

        elif (group == "obs"):
            cache_state(session_ID, key="# cells/obs", val=len(adata.index))
        elif (group == "var"):
            cache_state(session_ID, key="# genes/var", val=len(adata.index))
        with lock:
            store_store = zarr.DirectoryStore(zarr_cache_dir)
            store = zarr.open_group(store=store_store, mode='a')
            if (group in attribute_groups): # then -> write only that part of the object (fast)
                if (group == "var"):
                    if (np.nan in adata.var.index):
                        adata.var.index = pd.Series(adata.var.index).replace(np.nan, 'nanchung')
                        adata.var["gene_ID"] = pd.Series(adata.var["gene_ID"]).replace(np.nan, 'nanchung')
                        adata.var["gene_ids"] = pd.Series(adata.var["gene_ids"]).replace(np.nan, 'nanchung')
                # Write the requested attribute by reading any existing AnnData
                # and updating the attribute, then writing the full AnnData back
                # to the zarr store. This avoids reliance on anndata internals.
                try:
                    existing = ad.read_zarr(zarr_cache_dir) if os.path.exists(zarr_cache_dir) else None
                except Exception as e:
                    print(f"[WARN] ad.read_zarr existing failed: {e}")
                    try:
                        store_store = zarr.DirectoryStore(zarr_cache_dir)
                        store = zarr.open_group(store=store_store, mode='r')
                        try:
                            keys = list(store.group_keys())
                        except Exception:
                            try:
                                keys = list(store.keys())
                            except Exception:
                                keys = None
                        print(f"[DEBUG] zarr store keys: {keys}")
                    except Exception as e2:
                        print(f"[DEBUG] opening zarr store failed: {e2}")
                    existing = None

                if existing is None:
                    existing = ad.AnnData()

                # Try to set only the requested attribute; if that fails,
                # fall back to writing the full AnnData to avoid index/shape errors.
                try:
                    if group == 'obs':
                        existing.obs = adata if isinstance(adata, pd.DataFrame) else adata.obs
                    elif group == 'var':
                        existing.var = adata if isinstance(adata, pd.DataFrame) else adata.var
                    elif group == 'X':
                        existing.X = adata
                    elif group == 'layers':
                        existing.layers = adata
                    elif group == 'obsm':
                        existing.obsm = adata
                    elif group == 'varm':
                        existing.varm = adata
                    elif group == 'obsp':
                        existing.obsp = adata
                    elif group == 'varp':
                        existing.varp = adata
                    elif group == 'uns':
                        existing.uns = adata
                    elif group == 'raw':
                        existing.raw = adata
                    safe_write_zarr(existing, zarr_cache_dir)
                except Exception as e:
                    print(f"[WARN] group-specific zarr write failed: {e}")
                    try:
                        if isinstance(adata, ad.AnnData):
                            safe_write_zarr(adata, zarr_cache_dir)
                        else:
                            fallback = ad.AnnData()
                            try:
                                if group == 'obs':
                                    fallback.obs = adata if isinstance(adata, pd.DataFrame) else adata.obs
                                elif group == 'var':
                                    fallback.var = adata if isinstance(adata, pd.DataFrame) else adata.var
                                elif group == 'X':
                                    fallback.X = adata
                                elif group == 'layers':
                                    fallback.layers = adata
                            except Exception:
                                pass
                            try:
                                safe_write_zarr(fallback, zarr_cache_dir)
                            except Exception as e2:
                                print(f"[ERROR] writing zarr fallback also failed: {e2}")
                    except Exception as e3:
                        print(f"[WARN] inner fallback write block failed: {e3}")

                # write dense copies of X or layers if they're what was passed
                if (group == "X"):
                    dense_name = "X_dense"
                    write_dense.delay(zarr_cache_dir, "X",
                                      dense_name, chunk_factors)

                if (group == "layers"):
                    for l in list(adata.keys()): #layers was passed with parameter name "adata"
                        dense_name = "layers_dense/" + str(l)
                        write_dense.delay(zarr_cache_dir, "layers/" + l, 
                                          dense_name, chunk_factors)
                lock.release()
            else:
                # check that necessary fields are present in adata object
                if not ("leiden_n" in adata.obs):
                    if ("leiden" in adata.obs):
                        adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
                if not ("cell_ID" in adata.obs):
                    adata.obs["cell_ID"] = adata.obs.index
                if not ("cell_numeric_index" in adata.obs):
                    adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(0,len(adata.obs.index))))
                for i in ["user_" + str(j) for j in range(0, 6)]:
                    if not (i in adata.obs.columns):
                        adata.obs[i] = ["0" for k in adata.obs.index.to_list()]
                if not ("gene_ID" in adata.var):
                    adata.var["gene_ID"] = adata.var.index

                # make sure that there are no "nan" genes in the var index
                if (np.nan in adata.var.index):
                    adata.var.index = pd.Series(adata.var.index).replace(np.nan, 'nanchung')
                    adata.var["gene_ID"] = pd.Series(adata.var["gene_ID"]).replace(np.nan, 'nanchung')
                    adata.var["gene_ids"] = pd.Series(adata.var["gene_ids"]).replace(np.nan, 'nanchung')

                # Save the whole AnnData using public API. Construct or update
                # an existing AnnData when needed and write that back to zarr.
                try:
                    existing = ad.read_zarr(zarr_cache_dir) if os.path.exists(zarr_cache_dir) else None
                except Exception:
                    existing = None

                # Prefer writing the provided full AnnData directly.
                try:
                    if isinstance(adata, ad.AnnData):
                        safe_write_zarr(adata, zarr_cache_dir)
                    else:
                        existing = None
                        try:
                            existing = ad.read_zarr(zarr_cache_dir) if os.path.exists(zarr_cache_dir) else None
                        except Exception:
                            existing = None
                        if existing is None:
                            existing = ad.AnnData()
                        # assign attributes where available; guard each assignment
                        for attr in ['obs','var','obsm','varm','obsp','varp','uns','raw','X','layers']:
                            if hasattr(adata, attr):
                                try:
                                    setattr(existing, attr, getattr(adata, attr))
                                except Exception:
                                    pass
                        try:
                            safe_write_zarr(existing, zarr_cache_dir)
                        except Exception as e:
                            print(f"[ERROR] writing zarr fallback failed: {e}")
                except Exception as e:
                    print(f"[ERROR] writing adata.zarr directly failed: {e}")

                # making dense copies of X and layers (compressed to save disk space)
                dense_name = "X_dense"
                write_dense.delay(zarr_cache_dir, "X",
                                  dense_name, chunk_factors)

                for l in list(adata.layers.keys()):
                    dense_name = "layers_dense/" + str(l)
                    write_dense.delay(zarr_cache_dir, "layers/" + l, 
                                      dense_name, chunk_factors)

                lock.release()
            # set the file mod and access times to current time
            # then return adata as usual 
            os.utime(zarr_cache_dir) 
            return adata

def adata_cache_exists(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "adata_cache"
    zarr_cache_dir = filename  + ".zarr"
        
    if (os.path.exists(zarr_cache_dir) is True):
        return True
        
    return False

def adata_cache_group_exists(session_ID, group, store=None):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "adata_cache"
    zarr_cache_dir = filename  + ".zarr"
    
    if (store is None):
        try:
            store_store = zarr.DirectoryStore(zarr_cache_dir)
            store = zarr.open_group(store=store_store, mode='r')
            keys = list(store.group_keys())
            #store_store.close()
        except:
            return False
    else:
        keys = list(store.group_keys())
    
    if (group in keys):
        return True
    else:
        return False

def cache_gene_list(session_ID, gene_list=None):
    filename = save_analysis_path + str(session_ID) + "/gene_list_cache.pickle"
    lock_filename = filename + ".lock"
    lock = FileLock(lock_filename, timeout=20)
    
    if (gene_list is None):
        if (os.path.isfile(filename) is True):
            with open(filename, "rb") as f:
                gene_list = pickle.load(f)
        else:
            print("[ERROR] gene list cache does not exist at: " + str(filename))
            gene_list = None
        return gene_list
    else:
        gene_list.sort(key=str.lower)
        with lock:
            with open(filename, "wb") as f:
                pickle.dump(gene_list, f)
            return gene_list

# returns a list of cell_IDs 
# expects a list of lists of datapoint dictionaries
def get_cell_intersection(session_ID, list_of_selections,
                          pt_min=0, pt_max=1):
    obs = cache_adata(session_ID, group="obs")
    cell_intersection = set(obs.index.to_list())

    for cell_list in list_of_selections:
        if (cell_list in ["", 0, None, []]):
            continue
        
        cell_set = set()
        for cell in cell_list["points"]:
            cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
            cell_set.add(cell_ID)
        cell_intersection &= cell_set

    return cell_intersection

# returns dictionary of points that are in all violin selections
# across all genes that could be selected on
def get_violin_intersection(session_ID, violin_selected):
    if (violin_selected is None):
        return None

    # test which traces (genes) cells were selected from
    curves = set()
    for cell in violin_selected["points"]:
        curves.add(cell["curveNumber"])

    if (len(curves) == 0):
        # no cells selected - return None
        return None

    # get cells selected in each of these curves
    cells_in_curves = [set() for curve in curves]
    
    if (len(cells_in_curves) == 1):
        return violin_selected

    for cell in violin_selected["points"]:
        n = cell["curveNumber"]
        if (n in curves):
            cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
            (cells_in_curves[n]).add(cell_ID)

    # do the intersection
    cell_intersection = cells_in_curves[0]
    for cell_set in cells_in_curves:
        cell_intersection &= cell_set 

    # get the final list of points in dict format
    points = []
    for cell in violin_selected["points"]:
        cell_ID = (cell["text"]).rsplit(" ", 1)[-1]
        if (cell_ID in cell_intersection):
            points.append(cell)

    return {"points": points}

def get_ortholog_data(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "dmel_human_orthologs_disease_fb_2019_05.csv")

    ortholog_data = pd.read_csv(filename, sep="\t")
    ret = ortholog_data.loc[ortholog_data["Dmel_gene_symbol"] == selected_gene]
    return ret

def get_gene_snapshot(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "gene_snapshots_fb_2019_05.csv")

    snapshots = pd.read_csv(filename, sep="\t")
    ret = snapshots.loc[snapshots["GeneSymbol"] == selected_gene]
    #ret = ret["gene_snapshot_text"]
    return ret

def get_disease_data(session_ID, selected_gene):
    filename = (save_analysis_path 
    + "disease_model_annotations_fb_2019_05.csv")

    diseases = pd.read_csv(filename, sep="\t")
    ret = diseases.loc[diseases["Dmel_gene_ID"] == selected_gene]
    #ret = ret["gene_snapshot_text"]
    return ret

def cache_multicolor_scale(session_ID, multi_color_scale=None):
    if not (session_ID == None):
        filename = save_analysis_path + str(session_ID) + "/multi_color_scale.pickle"
    else:
        filename = save_analysis_path + "multi_color_scale.pickle"

    if (multi_color_scale is None):
        if (os.path.isfile(filename) is True):
            with open(filename, "rb") as f:
                multi_color_scale = pickle.load(f)
        else:
            print("[ERROR] multi_color_scale cache does not exist at: " + str(filename))
            multi_color_scale = MultiColorScale()
            with open(filename, "wb") as f:
                pickle.dump(multi_color_scale, f)
        return multi_color_scale
    else:
        with open(filename, "wb") as f:
            pickle.dump(multi_color_scale, f)
        return multi_color_scale

def remove_old_cache(n_days=1.5):
    n_sec_in_day = 86400
    max_time_in_sec = int(n_sec_in_day * n_days)
    
    now = time.time()

    for r,d,f in os.walk(save_analysis_path):
        for directory in d:
            timestamp = os.path.getmtime(os.path.join(r,directory))
            if (now-max_time_in_sec > timestamp):
                try:
                    print("[CLEANUP] removing " + str(os.path.join(r,directory)))
                    shutil.rmtree(os.path.join(r,directory))
                except Exception as e:
                    print("[ERROR] " + str(e))
                    pass

# implementation specific to custom zarr file cache with X_dense and layers_dense
def get_obs_vector(session_ID, var, layer="X"):
    save_dir = save_analysis_path  + str(session_ID) + "/"

    if (use_zarr is True):
        zarr_cache_dir = save_dir + "adata_cache" + ".zarr"
        if (os.path.exists(zarr_cache_dir) is True):
            try:
                adata = ad.read_zarr(zarr_cache_dir)
            except Exception as e:
                print(f"[ERROR] get_obs_vector: failed to read zarr: {e}")
                return None

            # prefer obs column if present
            if var in adata.obs.columns:
                return list(adata.obs[var])

            # try to find gene index by gene_ID column then var_names
            idx = None
            try:
                if 'gene_ID' in adata.var.columns:
                    idx = list(adata.var['gene_ID']).index(var)
            except Exception:
                idx = None

            if idx is None:
                try:
                    idx = list(adata.var_names).index(var)
                except Exception:
                    idx = None

            if idx is None:
                return None

            # extract vector from X or layers
            try:
                if layer == 'X':
                    vec = adata.X[:, idx]
                else:
                    vec = adata.layers[layer][:, idx]
                # convert to list
                try:
                    return list(vec)
                except Exception:
                    # If sparse matrix
                    return list(np.asarray(vec).ravel())
            except Exception as e:
                print(f"[ERROR] get_obs_vector: failed to extract vector: {e}")
                return None

def generate_marker_gene_table(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    uns = cache_adata(session_ID, group="uns")
    marker_genes = pd.DataFrame.from_records(uns["rank_genes_groups"]["names"])
    marker_genes.to_csv(save_dir + "marker_genes.csv")

def marker_genes_table_exists(session_ID):
    save_dir = save_analysis_path  + str(session_ID) + "/"
    filename = save_dir + "marker_genes.csv"
    if (os.path.isfile(filename) is True):
        return True
    return False

def to_rgba_string(rgb_tuple, opacity=1):
    ret = "rgba("
    for c in rgb_tuple:
        ret += str(c) + ","
    ret += str(opacity) + ")"
    return ret