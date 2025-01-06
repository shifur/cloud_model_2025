=== OVERVIEW ===

"MAIN.ipynb" is a boilerplate notebook to run ensemble simulations with the idealized cloud model of Posselt and Vukicevic (2010). It works by first sampling the model parameters from a uniform distribution (n_ens) and applying the model equations to each set of parameters. The parmap library allows this to be done in parallel (by running multiple processes at the same time).


=== COMPILATION ===

Before running the notebook, you will need to 1) install all necessary Python packages appearing in "MAIN.ipynb" and 2) compile the cloud_column_model code. 

For 2), first delete "cloud_column_model.x" from the top-level directory (it is a soft link) and inside the "cloud_column_model" directory (the actual exectutable). 

Then, following the instructions in the "README_derek.md" file, do "bash ./link_files.sh" and "bash ./make.sh 4". If compilation is successful, you should get the executable "cloud_column_model.x". Note that after running "bash ./link_files.sh", you might get the errors "ln:  ./dimensions.h: File exists", but you can safely neglect them.

After compiling, go back to the top-level directory and create a soft link to the cloud model executable (run "ln -sf cloud_column_model/cloud_column_model.x cloud_column_model.x").

You should be all set now and can run the main script "MAIN.ipynb".