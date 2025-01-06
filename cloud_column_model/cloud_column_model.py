"""
cloud_column_model.py -- a simple wrapper that encompasses the column cloud model first 
                         introduced in Posselt and Vukicevic (2010, MWR), and used in 
                         many subsequent microphysics parameter uncertainty and data 
                         assimilation experiments. 
                         
                         The model is written entirely in fortran, and reads a namelist 
                         (that should not need to be modified at all by this code). 
                         
                         The executable name is cloud_column_model.x, and it is expected that 
                         the executable will be in the current directory

                         There are three optional command line inputs to the model: 
                         1. The name of a file containing 11 microphysics parameters. 
                         This is a text file containing one row of 11 values, all of which 
                         are floating point.
                         2. The name of a text file containing a vector of output that is 36 elements long.
                         The vector is contained on a single row, and all values are floating point.
                         3. The name of the namelist file 
                         
                         This routine needs as input the vector of parameters that will be
                         input to the CRM, as well as the input file, output file, namelist file, and run number

                         Optional inputs include:
                         verbose (def False): Whether to write additional diagnostics to screen
                         params (def None): vector of parameter values. If provided, these will be written to the input file
                         rm_tmp (def True): whether to remove the temporary input and output files
"""

import subprocess
import uuid
# from parosse.framework.exewrap import ExeWrap
import numpy as np

class CRM1DWrap():
    def __init__(self, input_file, output_file, namelist_file, run_num=None, verbose=False, params=None, rm_tmp=True):
        """
        run_num: which run number this is - used in generating input and output files
        params: list of microphysical parameter values to write to the input file
        input_file: name and full path to input file
        output_file: name and full path to output file with run number appended
        rm_tmp: logical flag - whether to remove the input and output files
        verbose: logical flag - whether to print additional diagnostics to screen
        """
        self.infile_name = input_file
        if run_num is None:
            self.outfile_name = f'{output_file}-CRM1D-{uuid.uuid4().hex[:10]}'
        else:
            self.outfile_name = f'{output_file}_{run_num}'
        self.namefile_name = namelist_file
        self.run_num = run_num
        self.params  = params
        self.rm_tmp = rm_tmp
        self.verbose = verbose

    def pre(self):
        """
        Create input file for a cloud_column_model run
        """

        # Create temp input file
        pre_status_flag = False
        if self.params != None:
            if self.run_num is None:
                self.infile_name = f'{self.infile_name}-CRM1D-{uuid.uuid4().hex[:10]}'
            else:
                self.infile_name = f'{self.infile_name}_{self.run_num}'
            # Write parameter values into text file
            # First, convert list to string
            params_out_str = " ".join(str(item) for item in self.params)
            if self.verbose:
                print('Parameter values to be used in cloud_column_model.py: ',params_out_str)
            # Now, write to file using the "with ... as:" recommended method, which automagically closes the file
            with open(self.infile_name, 'w') as f_out:
              f_out.write(params_out_str)
            pre_status_flag = True
        # Return inputs
        return pre_status_flag

    def exe(self, inputs):
        """
        Do single run of cloud_column_model
        """
        pre_status_flag = inputs
        exe_status_flag = False

        if self.verbose:
            print('Input file name:    ',self.infile_name)
            print('Output file name:   ',self.outfile_name)
            print('Namelist file name: ',self.namefile_name)

        try:
            subprocess.check_call(['./cloud_column_model.x', self.infile_name, self.outfile_name, self.namefile_name])
            exe_status_flag = True
        except:
            print('CRM run failed!')
            exe_status_flag = False

        return exe_status_flag

    def post(self, inputs, outputs):
        """
        Read model output from file and delete temp files - input and output
        """
        pre_status_flag  = inputs
        exe_status_flag  = outputs
        post_status_flag = False

        # Read model output into list
        model_output_in = np.loadtxt(self.outfile_name)
        self.model_output = np.ndarray.tolist(model_output_in)

        # Remove temp files
        if (self.rm_tmp) & (exe_status_flag == True):
            subprocess.check_call(['rm', self.infile_name])
            subprocess.check_call(['rm', self.outfile_name])

        # If the pre-processor and exe have finished successfully, set post_status to true
        if pre_status_flag & exe_status_flag:
          post_status_flag = True

        return post_status_flag

    # def run(self, tagging=True):
    def __call__(self):
        """Default implementation for running the protocol steps:  pre, exe, post.
        """
        pre_status_flag = self.pre()
        exe_status_flag = self.exe(pre_status_flag)
        post_status_flag = self.post(pre_status_flag, exe_status_flag)

        return self.model_output, post_status_flag
