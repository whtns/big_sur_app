from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

import zarr
from numcodecs import Blosc
import scipy.sparse as sp

from .celery import task_queue
from filelock import FileLock
from progress_tracking import progress_tracker
from layouts import demo
if (demo == False):
    try:
        from configmodule.production_config import ProductionConfig as FlaskConfig
    except:
        from configmodule.default_config import DefaultConfig as FlaskConfig
else:
    from configmodule.default_config import DefaultConfig as FlaskConfig

@task_queue.task
def send_flask_mail(subject=None, sender=None,
                    recipients=None, body=None,
                    html=None):
    print("[DEBUG] running send_flask_mail")
    # Get mail config
    c = FlaskConfig()

    # Construct the message
    message = Mail(
        from_email=sender,
        to_emails=recipients,
        subject=subject,
        plain_text_content=body,
        html_content=html
    )
    try:
        sg = SendGridAPIClient(c.SENDGRID_API_KEY)
        print("[DEBUG] apikey: " + str(c.SENDGRID_API_KEY))
        response = sg.send(message)
        print(response.status_code)
        print(response.body)
        print(response.headers)
    except Exception as e:
        print(e)    
    return None

## Helper zarr caching
@task_queue.task
def write_dense(zarr_cache_dir, key, dense_name, chunk_factors):
    compressor = Blosc(cname='blosclz', clevel=3, shuffle=Blosc.SHUFFLE)
    lock_filename = str(zarr_cache_dir) + ".lock"
    lock = FileLock(lock_filename, timeout=60)
    with lock:
        store = zarr.open(zarr_cache_dir, mode='a')
        if (len(store[key]) == 3):
            # assume csr sparse matrix - parse as such
            array_keys = list(store[key].array_keys())
            X = sp.csr_matrix((store[key + "/" + array_keys[0]], 
                              store[key + "/" + array_keys[1]], store[key + "/" + array_keys[2]]))
        else:
            # assume dense matrix
            # TODO: checking for other cases of sparse matrices/mixed groups
            X = store[key]

        if ((not (dense_name in store))
        or  (X.shape != store[dense_name].shape)) :
            store.create_dataset(dense_name, shape=X.shape,
                                 dtype=X.dtype, fill_value=0, 
                                 chunks=(int(X.shape[0]/chunk_factors[0]), int(X.shape[1]/chunk_factors[1])),
                                 compressor=compressor, overwrite=True)
        if(sp.issparse(X) is True):
            X = X.tocoo()
            store[dense_name].set_coordinate_selection((X.row, X.col), X.data)
        else:
            store[dense_name] = X
    return None
### end celery queue task function definitions


## Progress Tracking Tasks

@task_queue.task
def track_progress(task_id, user_email, task_name, progress, total, details=None):
    """
    Track and store task progress, optionally sending email updates.
    
    Args:
        task_id: Unique identifier for the task
        user_email: Email of user running the task
        task_name: Name of the task being tracked
        progress: Current progress count
        total: Total items to process
        details: Optional details about current step
    """
    progress_tracker.update_progress(task_id, user_email, task_name, progress, total, details)
    return {'task_id': task_id, 'progress': progress, 'total': total}


@task_queue.task
def send_progress_notification(recipient_email, task_name, progress, total, 
                               details=None, task_status='in_progress'):
    """
    Send progress update email via SMTP.
    
    Args:
        recipient_email: Email address to send to
        task_name: Name of the task
        progress: Current progress count
        total: Total items to process
        details: Optional details about current step
        task_status: Status of task ('in_progress', 'completed', 'error')
    """
    return progress_tracker.send_progress_email(
        recipient_email, task_name, progress, total, details, task_status
    )


@task_queue.task
def send_task_completion(recipient_email, task_name, success=True, result_details=None):
    """
    Send task completion notification via SMTP.
    
    Args:
        recipient_email: Email address to send to
        task_name: Name of the completed task
        success: Whether task completed successfully
        result_details: Optional details about the result
    """
    return progress_tracker.send_completion_email(
        recipient_email, task_name, success, result_details
    )
### end progress tracking task definitions


# Correlation analysis tasks
@task_queue.task(bind=True)
def compute_correlations_task(self, session_ID, pvalue_threshold=0.001, correlation_threshold=0.1):
    """
    Compute gene correlations from HVG-selected genes in current session.
    Reports progress via task state.
    """
    import os
    from scipy.sparse import load_npz
    from helper_functions import cache_adata, save_analysis_path
    
    # Prefer local implementation first, fall back to external BigSuR package
    try:
        from correlations.correlation_functions import calculate_correlations
    except Exception:
        try:
            from BigSur.correlations import calculate_correlations
        except Exception:
            raise ImportError(
                "calculate_correlations is not available. Install the BigSuR package "
                "or ensure the local `src/correlations` module is present."
            )
    
    # Update progress: loading data
    self.update_state(state='PROGRESS', meta={'status': 'Loading session data...'})
    
    adata = cache_adata(session_ID)
    if adata is None:
        raise ValueError('No dataset in current session')
    
    # Check for HVG genes
    if 'highly_variable_user' not in adata.var.columns:
        raise ValueError('No HVG analysis found')
    
    hvg_mask = adata.var['highly_variable_user']
    if hvg_mask.sum() == 0:
        raise ValueError('No highly variable genes selected')
    
    # Update progress: subsetting data
    self.update_state(state='PROGRESS', meta={'status': f'Subsetting to {hvg_mask.sum()} HVGs...'})
    adata_hvg = adata[:, hvg_mask].copy()
    
    # Create output directory for correlation results
    session_dir = os.path.join(save_analysis_path, str(session_ID))
    os.makedirs(session_dir, exist_ok=True)
    correlation_output_dir = os.path.join(session_dir, 'correlation_results')
    os.makedirs(correlation_output_dir, exist_ok=True)
    
    # Update progress: computing correlations
    self.update_state(state='PROGRESS', meta={'status': 'Computing gene correlations...'})
    
    # Calculate correlations - results will be saved to disk
    # Note: write_out needs to end with / for BigSuR compatibility
    calculate_correlations(
        adata_hvg, 
        layer='counts',
        write_out=correlation_output_dir + '/',
        verbose=1
    )
    
    # Load the computed results from disk
    mcPCCs_path = os.path.join(correlation_output_dir, 'mcPCCs.npz')
    pvalues_path = os.path.join(correlation_output_dir, 'BH_corrected_pvalues.npz')
    
    # Check if files were created
    if not os.path.exists(mcPCCs_path) or not os.path.exists(pvalues_path):
            # List what files actually exist for debugging
            import glob
            actual_files = glob.glob(os.path.join(correlation_output_dir, '*'))
            error_msg = f'Correlation calculation did not produce output files. Files in {correlation_output_dir}: {actual_files}'
            raise ValueError(error_msg)
    
    # Verify we can load the results
    try:
        mcPCCs = load_npz(mcPCCs_path)
        pvalues = load_npz(pvalues_path)
    except Exception as e:
        raise ValueError(f'Failed to load correlation results: {str(e)}')
    
    if mcPCCs is None or mcPCCs.shape[0] == 0:
        raise ValueError('No correlation matrix returned')
    
    # Update progress: analyzing communities
    self.update_state(state='PROGRESS', meta={'status': 'Detecting communities and analyzing...'})
    
    return {
        'session_ID': session_ID,
        'mcPCCs_path': mcPCCs_path,
        'pvalues_path': pvalues_path,
        'n_genes': mcPCCs.shape[0],
        'status': 'completed'
    }