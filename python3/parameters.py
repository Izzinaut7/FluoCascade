 
# this is the cheap way to edit parameters without entering variables in the main python code.
# made to be imported in the main.py 

############################################
#
#       File parameters
#
#############################################

# put data in a directory
# the directory name can be relative or absolute path somewhere on Your computer
# This directory must exist before starting the absorption.py
datadir= "../../data"


# following variables are used to construct the name for the data files, so we can find and analyze them later:
appname="abs_09"


# some data can be reused, especially the transmission of the different fluorescence lines through the sample_thickness
# since this takes long, we can reuse it by referencing the previous calculation here:
simappname="abs_08"

#############################################
#
# Calculation parameters
#
#############################################

# V is for Voxel
Vsize=2E-7               #  m - I'd like 10nm, but no computer is going ot make it

#to test the code in a fast run, use larger voxel size, but take care that it is not larger than the box
Vsize=2E-6

# these describe the dimensions of the calculated volume
# it makes sense that the box size is larger than the beam diameter
sample_thickness=3.0E-4    # 0.3 mm
beam_dia = 3.2E-5 # 20 mu
# this is actually setting the half of the size of the box: the actual box is from -box_Y_size to box_Y_size:
box_Y_size = 30.0E-6 # 50 mu
box_Z_size = 17.0E-6 # 35 mu

# Intensity of the primary beam is used to get the dose values in a realistic order of magnitudes.
# the funny numer here is taken from the experimental data
I0=20334375 # photon/s - flux, beam intensity

# the list of primary beam energies to perform the calculation for. 
# one can run multiple energies or make a different parameter file for each
# E_beam=[5000,8000,12000,18000,25000, 30000,35000]
E_beam=[8000]


#############################################
#
#       Material parameters
#
#############################################

elements=["Ca", "P", "O", "N", "C"]


# TODO should use the database for this!
# but then I have several lines and all the proportionality factors...
# here I just assume that there is no K-beta, everything goes into one channel
fluorescence_E={
    "C":277,
    "N":392,
    "O":524,
    "P":2010,
    "Ca":3690
}

fluorescence_yield_K={
    "C":2.570e-03,
    "N":4.350e-03,
    "O":6.91E-3,
    "P":6.42E-2,
    "Ca":0.1687
}

# fluorescence Yield for Ca for L1 and L3 lines is:
# Shell  Yield     
# K   1.687e-01 
# L1  1.963e-04 
# L3  2.850e-04 
#
# Energy of the L-lines is 0.3-0.4 keV
#
