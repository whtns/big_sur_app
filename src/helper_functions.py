import os
import shutil
import time
import pickle
from datetime import datetime

from filelock import Timeout, FileLock
import zarr
from numcodecs import Blosc

# Lazy-import heavy scientific libs to avoid initializing JIT/llvmlite
# in the master process (pre-fork). Use `get_scanpy()` / `get_anndata()`
# to import inside worker runtime.
import pandas as pd
import numpy as np
import scipy as sp
_sc = None
_ad = None

def get_scanpy():
    global _sc
    if _sc is None:
        import scanpy as sc
        _sc = sc
    return _sc

def get_anndata():
    global _ad
    if _ad is None:
        import anndata as ad
        _ad = ad
    return _ad


def safe_read_zarr(zarr_cache_dir):
    """Read an AnnData from a zarr store with best-effort cleanup of helper
    dense groups that older helpers may have left behind (e.g., X_dense,
    layers_dense). If the first read fails, attempt to remove those groups
    and retry. If still failing, remove the zarr store and raise.
    """
    ad = get_anndata()
    try:
        return ad.read_zarr(zarr_cache_dir)
    except Exception as e:
        try:
            # Try to remove helper dense groups and retry
            if os.path.exists(zarr_cache_dir):
                try:
                    store = zarr.open_group(store=zarr.DirectoryStore(zarr_cache_dir), mode='a')
                    for g in ["X_dense", "layers_dense"]:
                        if g in store:
                            try:
                                del store[g]
                            except Exception:
                                pass
                except Exception:
                    pass
        except Exception:
            pass

        try:
            return ad.read_zarr(zarr_cache_dir)
        except Exception as e2:
            try:
                shutil.rmtree(zarr_cache_dir, ignore_errors=True)
            except Exception:
                pass
            raise
# Prefer public anndata I/O. Older code used anndata._io.zarr internals
# (read_dataframe/read_attribute/write_attribute) which are not stable.
# Use the public API (`ad.read_zarr`, `AnnData.write_zarr`) instead.
_ANNDATA_ZARR_INTERNALS = False

from status.status_functions import *
from plotting.multi_color_scale import MultiColorScale

from tasks.tasks import write_dense

# Paths and defaults for cached data
# Allow overriding via environment variables (container-friendly defaults)
save_analysis_path = os.environ.get('SAVE_ANALYSIS_PATH', '/app/data/cache/')
selected_datasets_path = os.environ.get('SELECTED_DATASETS_PATH', '/app/data/selected_datasets/')
user_dataset_path = os.environ.get('USER_DATASET_PATH', '/app/data/user_datasets/')

lock_timeout = 60
use_zarr = True

def generate_adata_from_10X(session_ID, data_type="10X_mtx"):
    """Load 10X-formatted data from the session raw_data directory."""
    data_dir = os.path.join(save_analysis_path, str(session_ID), "raw_data/")
    sc = get_scanpy()
    if data_type == "10X_mtx":
        adata = sc.read_10x_mtx(data_dir, cache=False)
    elif data_type == "10X_h5":
        adata = sc.read_10x_h5(os.path.join(data_dir, "data.h5ad"))
    else:
        print("[ERROR] data type not recognized - returning None")
        return None

    cache_adata(session_ID, adata)
    return adata

def load_selected_dataset(session_ID, dataset_key):
    """Load bundled sample dataset keyed by dropdown value."""
    dataset_dict = {
        "00000": "pbmc3k_raw",
        "00001": "pbmc3k_processed",
    }
    filename = dataset_dict.get(dataset_key)
    if filename is None:
        return None
    # Try several candidate locations so the app works both on-host and in-container.
    candidates = [
        os.path.join(selected_datasets_path, filename + ".h5ad"),
        os.path.join('/app/data', filename + ".h5ad"),
        os.path.join(save_analysis_path, filename + ".h5ad"),
        os.path.join('/Library/WebServer/Documents/BigSuR/selected_datasets', filename + ".h5ad"),
    ]

    path = None
    for p in candidates:
        try:
            if p and os.path.exists(p):
                path = p
                break
        except Exception:
            continue

    if path is None:
        print(f"[ERROR] selected dataset not found: {filename} (checked {candidates})")
        return None

    sc = get_scanpy()
    adata = sc.read_h5ad(path)

    state = {
        "filename": filename,
        "# cells/obs": len(adata.obs.index),
        "# genes/var": len(adata.var.index),
        "# counts": int(np.sum(adata.X)),
    }
    cache_state(session_ID, state)

    adata = cache_adata(session_ID, adata)
    return adata

def safe_write_zarr(adata_obj, zarr_cache_dir):
    """Safely write AnnData to a zarr store using atomic replace.

    - Write to a temp directory, then rename over the target.
    - Best-effort cleanup of helper groups on failure.
    """
    lock_filename = str(zarr_cache_dir) + ".lock"
    lock = FileLock(lock_filename, timeout=lock_timeout)

    ts = int(time.time() * 1000)
    pid = os.getpid()
    tmp_dir = f"{zarr_cache_dir}.tmp.{pid}.{ts}"
    backup_dir = f"{zarr_cache_dir}.bak.{pid}.{ts}"

    try:
        with lock:
            # clean any previous tmp/backup
            try:
                if os.path.exists(tmp_dir):
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                if os.path.exists(backup_dir):
                    shutil.rmtree(backup_dir, ignore_errors=True)
            except Exception:
                pass

            # write to tmp
            adata_obj.write_zarr(tmp_dir)

            # backup old and replace
            if os.path.exists(zarr_cache_dir):
                try:
                    os.rename(zarr_cache_dir, backup_dir)
                except Exception:
                    shutil.rmtree(zarr_cache_dir, ignore_errors=True)
            os.rename(tmp_dir, zarr_cache_dir)

            # cleanup backup
            try:
                if os.path.exists(backup_dir):
                    shutil.rmtree(backup_dir, ignore_errors=True)
            except Exception:
                pass
    except Exception as e:
        print(f"[WARN] safe_write_zarr failed: {e}")
        # fallback: remove helper groups then retry once
        try:
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)
            store = zarr.open_group(zarr_cache_dir, mode='a') if os.path.exists(zarr_cache_dir) else None
            if store:
                for g in ["X_dense", "layers_dense"]:
                    if g in store:
                        try:
                            del store[g]
                        except Exception:
                            pass
            # retry write
            adata_obj.write_zarr(zarr_cache_dir, overwrite=True)
        except Exception as e2:
            print(f"[ERROR] safe_write_zarr retry failed: {e2}")
            try:
                if os.path.exists(tmp_dir):
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                if os.path.exists(backup_dir):
                    shutil.rmtree(backup_dir, ignore_errors=True)
            except Exception:
                pass
            raise

def cache_adata(session_ID, adata=None, group=None, store_dir=None, store_name=None):
    """Read/write AnnData to the session zarr cache.

    If `adata` is None, read from cache; otherwise write to cache.
    `group` can be one of obs/var/obsm/varm/obsp/varp/layers/X/uns/raw for partial writes.
    """
    save_dir = store_dir if store_dir else os.path.join(save_analysis_path, str(session_ID), "")
    filename = (store_dir + store_name) if (store_dir and store_name) else os.path.join(save_dir, "adata_cache")
    chunk_factors = [3, 3] if store_dir else [150, 3]

    os.makedirs(save_dir, exist_ok=True)

    zarr_cache_dir = filename + ".zarr"
    attribute_groups = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "uns", "raw"]

    lock_filename = os.path.join(save_dir, "adata.lock")
    lock = FileLock(lock_filename, timeout=lock_timeout)

    if adata is None:
        if not os.path.exists(zarr_cache_dir):
            return None
        def _clean_dense_groups(path):
            try:
                store = zarr.open_group(store=zarr.DirectoryStore(path), mode='a')
                for g in ["X_dense", "layers_dense"]:
                    if g in store:
                        try:
                            del store[g]
                        except Exception:
                            pass
            except Exception:
                pass

        try:
            adata = safe_read_zarr(zarr_cache_dir)
        except Exception as e:
            print(f"[ERROR] ad.read_zarr failed: {e}")
            # Clean up helper dense groups then retry once
            _clean_dense_groups(zarr_cache_dir)
            try:
                adata = safe_read_zarr(zarr_cache_dir)
            except Exception as e2:
                print(f"[ERROR] ad.read_zarr retry failed: {e2}")
                try:
                    shutil.rmtree(zarr_cache_dir, ignore_errors=True)
                    print(f"[WARN] deleted stale zarr cache at {zarr_cache_dir}")
                except Exception as e_rm:
                    print(f"[WARN] failed to delete stale zarr cache: {e_rm}")
                return None

        if group in attribute_groups:
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
        return adata

    # write path
    if group is None:
        cache_state(session_ID, key="# cells/obs", val=len(adata.obs.index))
        cache_state(session_ID, key="# genes/var", val=len(adata.var.index))
        if ("total_counts" in adata.obs):
            cache_state(session_ID, key="# counts", val=int(np.sum(adata.obs["total_counts"])))
        else:
            cache_state(session_ID, key="# counts", val=int(np.sum(adata.X)))
    elif group == "obs":
        cache_state(session_ID, key="# cells/obs", val=len(adata.index))
    elif group == "var":
        cache_state(session_ID, key="# genes/var", val=len(adata.index))

    with lock:
        if group in attribute_groups:
            try:
                existing = safe_read_zarr(zarr_cache_dir) if os.path.exists(zarr_cache_dir) else get_anndata().AnnData()
            except Exception:
                existing = get_anndata().AnnData()

            # assign group-specific data
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

            if group == "X":
                write_dense.delay(zarr_cache_dir, "X", "X_dense", chunk_factors)
            if group == "layers":
                for l in list(adata.keys()):
                    write_dense.delay(zarr_cache_dir, f"layers/{l}", f"layers_dense/{l}", chunk_factors)
        else:
            # ensure required fields exist
            if "leiden_n" not in adata.obs and "leiden" in adata.obs:
                adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])
            if "cell_ID" not in adata.obs:
                adata.obs["cell_ID"] = adata.obs.index
            if "cell_numeric_index" not in adata.obs:
                adata.obs["cell_numeric_index"] = pd.to_numeric(list(range(len(adata.obs.index))))
            for i in ["user_" + str(j) for j in range(0, 6)]:
                if i not in adata.obs.columns:
                    adata.obs[i] = ["0" for _ in adata.obs.index]
            if "gene_ID" not in adata.var:
                adata.var["gene_ID"] = adata.var.index

            safe_write_zarr(adata, zarr_cache_dir)

            write_dense.delay(zarr_cache_dir, "X", "X_dense", chunk_factors)
            for l in list(adata.layers.keys()):
                write_dense.delay(zarr_cache_dir, f"layers/{l}", f"layers_dense/{l}", chunk_factors)

        # touch the path to update mtime
        try:
            os.utime(zarr_cache_dir, None)
        except Exception:
            pass

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
                adata = safe_read_zarr(zarr_cache_dir)
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