#!/bin/bash

# Runs "make" in each directory necessary to run the column model
# Takes as input a numerical value, which is used to determine 
# which directories to compile in. Default is all, but user can specify:
#
# 1: main directory only
# 2: main and microphysics
# 3: main, microphysics, and dynamics
# 4: main, microphysics, dynamics, and radiation (default)

let ctype=${1}

# Microphysics
if [ $ctype -ge 2 ]; then
echo ''
echo '***************************************'
echo '        Compiling microphysics'
echo '***************************************'
echo ''
# rm -f microphysics.a
cd microphysics_29Sept2008
# # cd microphysics_16July2010
# # cd microphysics_03Aug2010
# cd microphysics_Steward_Feb2013
make clean; make
cp -p microphysics.a ../
cd ../
# ln -s microphysics_03Aug2010/microphysics.a .
# ln -s microphysics_29Sept2008/microphysics.a .
fi

# Dynamics
if [ $ctype -ge 3 ]; then
echo ''
echo '***************************************'
echo '        Compiling dynamics'
echo '***************************************'
echo ''
# rm -f dynamics.a
cd dynamics_29Sept2008
# cd dynamics_Steward_Feb2013
make clean; make
cp -p dynamics.a ../
cd ../
# ln -s dynamics_29Sept2008/dynamics.a .
fi

# Radiation
if [ $ctype -ge 4 ]; then
echo ''
echo '***************************************'
echo '        Compiling radiation'
echo '***************************************'
echo ''
# rm -f radiation.a
cd radiation/
# # cd radiation_16July2010/
# # cd radiation_03Aug2010/
# cd radiation_Steward_Feb2013
make clean; make
cp -p radiation.a ../
cd ../
# ln -s radiation/radiation.a .
# ln -s radiation_03Aug2010/radiation.a .
# rm -f radiation_2009.a
cd radiation_2009/
make clean; make
cp -p radiation_2009.a ../
cd ../
# ln -s radiation_2009/radiation_2009.a .
fi

# Create main executable from libraries
if [ $ctype -ge 1 ]; then
echo ''
echo '***************************************'
echo '        Compiling main'
echo '***************************************'
echo ''
make clean; make
fi

# If called for, just clean
if [ $ctype -eq 0 ]; then

rm -f microphysics.a
cd microphysics_29Sept2008
make clean
cd ../

rm -f dynamics.a
cd dynamics_29Sept2008
make clean
cd ../

rm -f radiation.a
cd radiation/
make clean
cd ../

rm -f radiation_2009.a
cd radiation_2009/
make clean
cd ../

make clean

fi