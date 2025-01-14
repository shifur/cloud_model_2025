"""
This function creates an ensemble of variables for use in simulations, data assimilation, etc.

It takes as input:
  n_ens:        The total desired number of ensemble members
  p1, p2:       parameters of the prior distribution (e.g., mean and std dev of Normal)
  pmin, pmax:   the min/max of the parameter values
  ens_gen:      the method of ensemble generation (uniform, normal, LHS, Sobol, Halton, Hammersly)

It returns a python list (param_list) dimensioned (n_ens, n_params) 
that can then be used in parmap to drive an ensemble of model simulations

The generation scheme is set by "ens_gen", and available schemes include:

Samples drawn from a probability distribution, including:
  uniform:    samples drawn from within a bounded uniform distribution with min = pmin, max = pmax
  normal:     samples drawn from a Normal distribution with mean = p1 and std dev = p2
  gamma:      samples drawn from a gamma distribution with mean = p1 and std dev = p2. 
              k,theta are computed from mean and std dev as k = mu^2/sig^2, theta = mu^2/sig

Quasi-random samples generated by the scikit-optimize routines:
  LHS:        latin hypercube sampling
  Sobol:      Sobol sequence
  Halton:     Halton sequence
  Hammersly:  Hammersly sequence

Full grid, also generated by the scikit-optimize routine:
  Grid:       Grid search. Will generate an ensemble length that is N^N with N the largest root of n_ens
              So, if the user wants M bins for each of N parameters, they must set n_ens = M**N
              Note that there are options to "Grid": 
                border = 
                  "include" - Always include pmin and pmax as points in the list
                  "exclude" - Never include pmin and pmax as points in the list
                  "only" - Only include pmin and pmax as points in the list
                use_full_layout = 
                  When True, a  full factorial design is generated and
                  missing points are taken from the next larger full factorial
                  design, depending on `append_border`
                  When False, the next larger  full factorial design is
                  generated and points are randomly selected from it.

Documentation on these generators can be found at 

https://scikit-optimize.github.io/stable/

and the code is here:

https://github.com/scikit-optimize/scikit-optimize/blob/master/skopt/sampler/

and examples of how each algorithm performs can be found here:

https://scikit-optimize.github.io/stable/auto_examples/sampler/initial-sampling-method.html

Derek Posselt
JPL
29 March 2021

Changes:
17-Feb-2023, DJP: Added the option to generate an ensemble using Morris screening. 
                  Requires the SALib; https://salib.readthedocs.io/en/latest/index.html
                  For Morris to work, a dictionary of inputs is required, containing:
                      moat_dict = {
                                  'num_vars': nparam,
                                  'names': pnames,
                                  'bounds': pbounds
                                  }
                      names is a list of strings, each associated with one parameter
                      bounds is an N-D numpy array dimensioned (num_vars,2) and containing min and max value for each parameter
                      The user should also provide a number of Morris "levels" - defaults to 10, but can range from 4 to 20
                        The number of levels are the number of bins each parameter is divided into - it effectively sets the dx perturbation length
                          Larger = more efficient spanning of the space
                          Smaller = more fine grained analysis
                        It is not possible to set different levels for different parameters, but the fine-grained-ness of the analys
                          can be controlled by adjusting the parameter range via pmin and pmax
                      Note that the number of ensemble members will be larger than n_ens

"""

# Import necessary modules
# import time
import sys
import numpy as np

# Scikit-Optimize sampling methods
from skopt.space import Space
from skopt.sampler import Sobol
from skopt.sampler import Lhs
from skopt.sampler import Halton
from skopt.sampler import Hammersly
from skopt.sampler import Grid

def create_ensemble ( n_ens, p1, p2, pmin, pmax, ens_gen, pnames=None, num_levels=10 ):

    # n_ens = self.n_ens
    # p1 = self.p1
    # p2 = self.p2
    # pmin = self.pmin
    # pmax = self.pmax
    # ens_gen = self.ens_gen

    # If using Morris screening, try importing the module and then create the dictionary
    if ens_gen == 'Morris':
      try:
        import SALib.sample.morris as morris_sample
      except:
        print('User requested Morris sampling')
        print('Could not import SALib.sample.morris')
        print('Returning to calling function')
        return

      # Create the Morris dictionary
      pbounds = np.transpose(np.array([pmin,pmax]))
      nparam = len(pmin)
      # If the user has not provided a list of names, generate it
      if pnames is None:
        pnames = []
        for p in range(nparam):
          pnames.append('p'+str(p))
      
      moat_dict = {
                  'num_vars': nparam,
                  'names': pnames,
                  'bounds': pbounds
                  }

      param_array = morris_sample.sample(moat_dict, n_ens, num_levels)

      # Reset the values of n_ens
      n_ens = len(param_array[:,0])

      # Cast the parameter array into a list
      param_list = param_array.tolist()

    # If using any of the scikit-optimize routines, set up the Space with min and max values
    # and then go ahead and generate all of the parameters as one long list
    elif ens_gen == 'LHS' or ens_gen == 'Sobol' or ens_gen == 'Halton' or ens_gen == 'Hammersly' or ens_gen == 'Grid':

      space = Space(np.transpose(np.array([pmin,pmax])))

      # Set up the perturbation generator
      if ens_gen == 'LHS':
        pert_gen = Lhs(lhs_type="classic", criterion=None)
      elif ens_gen == 'Sobol':
        pert_gen = Sobol()
      elif ens_gen == 'Halton':
        pert_gen = Halton()
      elif ens_gen == 'Hammersly':
        pert_gen = Hammersly()
      elif ens_gen == 'Grid':
        pert_gen = Grid(border="include", use_full_layout=True) 

      # Every generator uses the same command, and returns a list
      param_list = pert_gen.generate(space.dimensions, n_ens)

    # Otherwise, random generation
    # Uniform
    elif ens_gen == 'uniform':

      param_list = np.random.uniform(pmin, pmax, (n_ens,len(pmin))).tolist()
      
    # Normal, constrained to lie within pmin,pmax...
    elif ens_gen == 'normal':
      param_list = np.empty((n_ens,len(pmin)))
#       for p in range(len(pmin)):
#         param_list[:,p] = np.random.normal(p1[p], p2[p], n_ens)
      for nn in range(n_ens):
        # Loop over parameters
        for p in range(len(pmin)):
          # Ensure parameter value is in range between pmin and pmax
          p_pert = np.random.normal(p1[p], p2[p]) # initialize
          # Find a suitable value
          while p_pert < pmin[p] or p_pert > pmax[p]:
            p_pert = np.random.normal(p1[p], p2[p])
          # Fill the parameter array
          param_list[nn,p] = p_pert
      # Convert to list
      param_list = param_list.tolist()

    # Gamma, could implement again constraint to lie within pmin,pmax, but not for now
    elif ens_gen == 'gamma':
      # Calculate k and theta for gamma from parameter mean (p1) and sigma (p2)
      kgam  = np.array(p1)[:]**2 / np.array(p2)[:]**2
      thgam = np.array(p2)[:]**2 / np.array(p1)[:]
      param_list = np.empty((n_ens,len(pmin)))
      # Loop over parameters
      for p in range(len(pmin)):
        # Create Gamma ensemble for this parameter
        param_list[:,p] = np.random.gamma(np.tile(kgam[p],n_ens),np.tile(thgam[p],n_ens))
      # If we want constraint...
#       for nn in range(n_ens):
#         # Loop over parameters
#         for p in range(len(pmin)):
#           # Ensure parameter value is in range between pmin and pmax
#           p_pert = pmin[p] - 1.0 # initialize
#           # Find a suitable value
#           while p_pert < pmin[p] or p_pert > pmax[p]:
#             p_pert = np.random.gamma(kgam[p], thgam[p])
#           # Fill the parameter array
#           param_list[nn,p] = p_pert
      # Convert to list
      param_list = param_list.tolist()    

    # If the user has requested a generation scheme we do not recognize, exit.
    else:
      print('User requested invalid prior ensemble generation type in create_ensemble.py: ',ens_gen)
      sys.exit("Stopping.")

    return param_list

# def __call__( self, n_ens, p1, p2, pmin, pmax, ens_gen ):

#   """
#   Main calling function
#   """
#   self.n_ens = n_ens
#   self.p1 = p1
#   self.p2 = p2
#   self.pmin = pmin
#   self.pmax = pmax
#   self.ens_gen = ens_gen

#   param_list = create_ensemble()

#   return param_list

