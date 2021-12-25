# FluoCascade
One-time written Pyhon calculation of the fuorescence cascades in a bone to estimate the bone damage in some special cases.

Purpose:

Calculate the secondary fluorescence in bone segment during the radiation with X-rays. For a given incident energy and beam profile it calculates the locally absorbed number of photons and dose. It uses volume elements (voxels) to keep track of local variables - hence large ammount of physical memory and HD space is required to run it.

Elements.hdf

Data file I collected from DABAX tables, mostly used for the scattering cross sections. The use of this is deprecated, and the next version is going to use xraylib module instead of this. (https://github.com/tschoonj/xraylib)

absorption.py

After calculating the number of absorbed photons, it calculates the number of photons absorbed by all elements in the voxel, as well as the total ammount of the fluorescence photons originating rom the voxel. 

Output:
the generated fields are exported into npz file format, the NumPy compressed arrays. To examine these arrays use plot.py and save_txy.py.

parameters.py is used to modify the calculation parameters, as well as the data directory. 

Note: save_txt does not work yet.
