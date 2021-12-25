# FluoCascade
One-time written Pyhon calculation of the fuorescence cascades in a bone to estimate the bone damage in some special cases.

Purpose:

Calculate the secondary fluorescence in bone segment during the radiation with X-rays. For a given incident energy and beam profile it calculates the locally absorbed number of photons and dose. It uses volume elements (voxels) to keep track of local variables - hence large ammount of physical memory and HD space is required to run it.

The code is not intended to be used on other problems, although I plan to do it in some later version. Specific properties of the material in question, namely dry bone are deeply hard-coded in the scripts. 

# Files inthe repository:

Elements.hdf
Data file I collected from DABAX tables, mostly used for the scattering cross sections. The use of this is deprecated, and the next version is going to use xraylib module instead of this. (https://github.com/tschoonj/xraylib)

absorption.py
After calculating the number of absorbed photons, it calculates the number of photons absorbed by all elements in the voxel, as well as the total ammount of the fluorescence photons originating rom the voxel. 

Output:
the generated fields are exported into npz file format, the NumPy compressed arrays. To examine these arrays use plot.py and save_txy.py.

Note: save_txt does not work yet.

# Input parameters: 
parameters.py is used to modify the calculation parameters, as well as the data directory. 

Important: data directory must exist before starting the simulation.

For further description read the comments in the parameters.py
Program is tested with 
Debian linux, Python 3.7
Windows10, Anaconda3, Python 3.8

# Installation Instructions

Download the code from the repository

Create a directory for the results and assign the location to the ddatadir variable in "parameters.py"

Edit other calcualtion parameters in the "parameters.py"

In a terminal (Anaconda Prompt in Windows) change to the directory where the python code is

Start the calculateion by typing
python3 absorption.py

(in Anaconda you should be in an active Python 3 environment, and be able to type "python absorption.py" without explicitely typing the version)






