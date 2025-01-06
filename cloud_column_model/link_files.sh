#!/bin/bash

# This script creates soft links to the include files necessary for code compilation in the sub directories.
#
# Necessary files are:
# dynamics_29Sept2008: dimensions.h
# microphysics_29Sept2008: dimensions.h
# radiation: dimensions.h, radiation.h
# radiation_2009: dimensions.h

# Set home directory
#homedir=`pwd` # hch change (5 March 2024)
homedir='../' # hch change (5 March 2024)

# Set include directory
incldir='../include/'

# Set file names
dims='dimensions.h'
rad='radiation.h'

# Create soft links
cd ./dynamics_29Sept2008/
ln -s ${incldir}${dims} .
cd ${homedir}/microphysics_29Sept2008/
ln -s ${incldir}${dims} .
cd ${homedir}/radiation/
ln -s ${incldir}${dims} .
ln -s ${incldir}${rad}  .
cd ${homedir}/radiation_2009/
ln -s ${incldir}${dims} .

# Go home
#cd ${homedir} # hch change (5 March 2024)

