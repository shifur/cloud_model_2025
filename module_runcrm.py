'''
Import this module in notebooks to enable the parmap capability.

---
Author: Hristo Chipilski, 29 Apr 2024.
'''

from cloud_column_model import cloud_column_model

def runcrm(runvar):
    """
    Run a single 1D CRM simulation.

    Input
    ----------
    runvar : tuple
        Parameters necessary for a cloud_colum_model run.

    Returns
    -------
    model_output: 
        Numpy vector of output variables from the CRM simulation
    """
    input_file, output_file, namelist_file, run_num, params = runvar
    crm1d = cloud_column_model.CRM1DWrap(input_file, output_file, namelist_file, run_num=run_num, params=params)
    model_output,crm_status = crm1d()
    return model_output
