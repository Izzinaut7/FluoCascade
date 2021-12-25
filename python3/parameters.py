 
# this is the cheap way to edit parameters without entering variables in the main python code.
# made to be imported in the main.py 

#####################
# file parameters:

datadir = "/data-zmk/basic-res-shared/tomo/Katrein/3rd_paper_DAMAGE/IVO_SIMULATION/ALLE_DATEN/secondary_emission"

# put data somewhere temporarilly
# must have this directory before starting the main.py
datadir= "../../data"

# following variables are used to construct the name for the data files, so we can find and analyze them later:
appname="abs_09"
#simappname1="in_src_001"   # loads the data from raytracing simulation
simappname1="sec_abs_001"  # loads the data from raytracing simulation, old version
simappname2="abs_08"      # loads the data we calculated in a previous run. If we didn't calculate anythiong, set this to appname


# V is for Voxel
Vsize=2E-7               #  m - I'd like 10nm, but no computer is going ot make it
#to test the code, use larger voxel size:
Vsize=2E-6

sample_thickness=3.0E-4    # 0.3 mm
beam_dia = 3.2E-5 # 20 mu
box_Y_size = 30.0E-6 # 50 mu
box_Z_size = 17.0E-6 # 35 mu

I0=20334375 # photon/s - flux, beam intensity


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

E_beam=[5000,8000,12000,18000,25000, 30000,35000]
E_beam=[5000,12000,25000]
E_beam=[3700]