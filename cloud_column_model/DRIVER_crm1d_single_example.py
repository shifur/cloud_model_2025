#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
crm1d_single_example.py -- This Python 3.x code simply runs a single experiment with the 1D column model
                           It does a few things that a full experiment would be expected to do:
                           1. Sets a "true" parameter list and runs the model to produce reference data
                           2. Generates a parameter set randomly from within a bounded uniform distribution
                           3. Runs the model with the new parameters
                           4. Prints output to screen
                           
                           Requires: cloud_column_model.py -- the wrapper for cloud_column_model.x

"""

# Import necessary modules, including the column model wrapper
import numpy as np
import cloud_column_model

# Set up input and output file names and initialize the run number
input_file = 'run_one_crm1d.txt'
output_file = 'crm1d_output.txt'
namelist_file = 'namelist_3h_t30-180.f90'
run_num = 0

# Set the random number seed for reproducibility
rand_seed = 42

# Set the bounds for draws of new parameter values
p_min = [  63.0,    0.10,     75.0,    0.1,  0.0001,   0.001,   0.001,    0.05,    0.05,    1.e-6,  1.e-6 ]
p_max = [1000.0,     1.0,   1200.0,   0.85,     5.0,     5.0,     5.0,     1.0,     1.0,    3.e-3,  2.e-3 ]

# Construct a reference list of parameter values and write to input file
p_ref = [200.0,  0.3,  400.0,  0.4,   0.5,    0.5,   0.5,   0.2,   0.4,  1.e-3, 6.e-4]
# Create a string from the list
params_out_str = " ".join(str(item) for item in p_ref)
# Now, write to file using the "with ... as:" recommended method, which automagically closes the file
with open(input_file, 'w') as f_out:
  f_out.write(params_out_str)

# Now, run the model
# crm1d = cloud_column_model.CRM1DWrap('ref_'+input_file, 'ref_'+output_file, namelist_file, run_num=run_num, params=p_ref, rm_tmp=False)
crm1d = cloud_column_model.CRM1DWrap('ref_'+input_file, 'ref_'+output_file, namelist_file, params=p_ref, rm_tmp=True)
# model_output_ref = crm1d.model_output
model_output_ref = crm1d()

# Print output
print('Reference input to model:    ',p_ref)
print('Reference output from model: ',model_output_ref)
print('')

# Now, generate parameters randomly from a bounded uniform distribution
# Set the seed
np.random.seed(rand_seed)

# Generate the parameters
nparam = len(p_ref)
p_pert_vec = np.random.uniform(p_min, p_max, nparam)
# Convert np vector to list
p_pert = np.ndarray.tolist(p_pert_vec)

# Now, run the model
crm1d = cloud_column_model.CRM1DWrap('pert_'+input_file, 'pert_'+output_file, namelist_file, run_num=run_num, params=p_pert, rm_tmp=False)
# crm1d.run(tagging=False)
# model_output_pert = crm1d.model_output
model_output_pert = crm1d()

# Print output
print('Perturbed parameters:        ',p_pert)
print('Output from model with pert: ',model_output_pert)

