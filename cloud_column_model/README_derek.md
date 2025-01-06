# cloud_column_model
This repository contains code for the cloud column (vertical 1D) model described in Posselt and Vukicevic (2010, Monthly Weather Review).

All code is written in Fortran (primarily fortran 77), and has been tested with the GNU fortran compiler (version 8.3) on a macbook pro. 

GNU fortran compilers can be obtained from the Mac OS X HPC sourceforge site at http://hpc.sourceforge.net/

All experiment settings can be found in the namelist file: namelist.f90

To run the code, first link the include files by issuing "bash ./link_files.sh"
Then, compile using "bash ./make.sh 4". The argument to the makefile determines which sub-directories are compiled; "4" compiles everything.
The executable is named cloud_column_model.x 

The default experiment runs a single simulation of the column model. The output to screen should be equivalent to what is provided in the file "test_results.txt"

Note that all units in the model are cgs.

## Model Structure

The top level driver (main routine) is cloud_column_model.f90. This reads the namelist, configures the experiment, and calls a routine to run an experiment. The default is to run a single model integration. The main model driver is drive_micro_rad.f90, which integrates the model over all timesteps, calling the dynamics (vertical motion only in this case), microphysics, and radiation. 

Dynamics code should not need any modification, and is all located in dynamics_29Sept2008/. 

The radiation code should also not need modification, and is all located in radiation/. There is a second radiation code "radiation_2009", but this is not currently supported. Note that radiation is run with a different frequency set in the namelist as "rad_dt". If this is set > model run time, then radiation will only be called when it is needed for "simulated observations" (e.g., outgoing longwave and shortwave radiation).

The microphysics code is located in microphysics_29Sept2008/. There are two key microphysics routines. consatrh.for, which sets up the microphysics, is run once at the beginning of the integration. saticerh.for contains the main microphysics code, and is run once per timestep. 

## Namelist

The namelist is comprised of 5 sections (records).
RECORD0: Model run parameters, including: run length; timestep; vertical grid spacing; and temperature, water vapor, and vertical velocity forcing.
RECORD1: Microphysics parameter minimum and maximum values
RECORD2: Settings to control which times are output, and which state variables are used as "observations"
RECORD3: Microphysics parameter reference (true) value and flag to determine whether each parameter is perturbed
RECORD4: Settings to control ensemble simulations

Note that the names of the parameters in the namelist correspond to the names in the microphysics code, not the publications. Definitions are as follows:
1. as = a_s, the coefficient in the snow fallspeed-diameter relationship
2. bs = b_s, the exponent in the snow fallspeed-diameter relationship
3. ag = a_g, the coefficient in the graupel fallspeed-diameter relationship
4. bg = b_g, the exponent in the graupel fallspeed-diameter relationship
5. tnw = N_0r, the intercept parameter in the (exponential) rain particle size distribution
6. tns = N_0s, the intercept parameter in the (exponential) snow particle size distribution
7. tng = N_0g, the intercept parameter in the (exponential) graupel particle size distribution
8. roqs = rho_s, the snow density
9. roqg = rho_g, the graupel density
10. bnd21 = q_c0, the cloud-to-rain autoconversion threshold
11. bnd1 = q_i0, the ice-to-snow autoconversion threshold (not used)

