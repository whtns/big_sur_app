import numpy as np
def check_csr_data_for_int_floats(csr_matrix_obj):
    """
    Checks if the non-zero data elements in a csr_matrix are floats
    that represent integer values.
    """
    if not np.issubdtype(csr_matrix_obj.data.dtype, np.floating):
        print("Data type is not float, skipping check.")
        return True # Or False, depending on desired behavior

    for value in csr_matrix_obj.data:
        if value != 0.0 and not value.is_integer():
            return False  # Found a non-zero float that is not an integer
    return True  # All non-zero floats represent integers