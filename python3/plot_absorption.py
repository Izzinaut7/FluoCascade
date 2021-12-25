#!/usr/bin/env  python3

__author__ = "Ivo Zizak (ivo.zizak@helmholtz-berlin.de)"
__version__ = "$ 0.09 $"
__date__ = "$Date: 2021/12/23 00:00:00 $"
__copyright__ = "Copyright (c) 2021 Ivo Zizak"
__license__ = "Python"


import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import ndimage   # for convolve
from scipy import signal   # for convolve
import pandas as pd
import re
# this is needed only for testing:
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import h5py
import time

from parameters import *

#appname="abs_08"

for E in E_beam:
            
    filename="{}/{}_F3D_{}_{:04d}.npz".format(datadir,appname,"I0",E)
    # files:
    #  absorption=F3D_I0_abs, 
    #  dose=dose, 
    #  dose_Col=dose_Col, 
    #  dose_HAP=dose_HAP )
    I0_z=np.load(filename)

    filename="{}/{}_F3D_{}_{:04d}.npz".format(datadir,appname,"Tot",E)
    Tot_z=np.load(filename)

    ### construct axes for plots
    dose=I0_z["dose"]
    absorption = I0_z["absorption"]
    dose_HAP = I0_z["dose_HAP"]
    dose_Col = I0_z["dose_Col"]

    Tot_dose=Tot_z["dose"]
    Tot_dose_Col=Tot_z["dose_Col"]

    # make axes
    Nx=dose.shape[0]
    Ny=dose.shape[1]
    Nz=dose.shape[2]
    Vsize=2E-7 * 1E6 # so we have everything in mu
    
    x_axis = np.linspace(0, Vsize*Nx, Nx)
    y_axis = np.linspace(-Vsize*(Ny-1)/2, Vsize*(Ny-1)/2, Ny)
    z_axis = np.linspace(-Vsize*(Nz-1)/2, Vsize*(Nz-1)/2, Nz)
    print(x_axis)
    xs,ys =np.meshgrid( y_axis, x_axis)

    
    print("I0 dose  . . . . . . {:g}".format(dose.sum()))
    print("Collagen I0 dose . . {:g}".format(dose_Col.sum()))
    print("HAP I0 dose  . . . . {:g}".format(dose_HAP.sum()))
    print("Total dose . . . . . {:g}".format(Tot_dose.sum()))
    print("Total Col dose . . . {:g}".format(Tot_dose_Col.sum()))
    print("Total Col dose (avg) {:g}".format(Tot_dose_Col.mean()))

    depth=40
    depthmu=depth*Vsize
    
    
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title    ##
    
    
    # in case you import matplotlib:
    #import matplotlib
    #SMALL_SIZE = 8
    #matplotlib.rc('font', size=SMALL_SIZE)
    
    abs_el={}
    for el in elements:
        filename="{}/{}_F3D_{}_{:04d}.npz".format(datadir,appname,el,E)
        abs_el[el]=np.load(filename)
        # absorption = F3D_Fluo_Abs, 
        # fluorescence = F3D_Fluo_elem, 
        # dose=dose_fluo, 
        # dose_col=dose_fluo_Col
        fig_el, axes = plt.subplots(2,2,figsize=(14,10), dpi=300)
        fig_el.suptitle("Abs. {} E={} depth={:.1f}".format(el,E, depthmu))
        
        axes[0,0].set_title("Intensity in ph/s / voxel")
        axes[0,0].plot(y_axis, abs_el[el]["fluorescence"][depth,:,int((Nz-1)/2)]  , label="fluorescence" )
        axes[0,0].plot(y_axis, abs_el[el]["absorption"][depth,:,int((Nz-1)/2)]  , label="absorption" )
        axes[0,0].legend()

        #axes[1,0].set_title("Absorption")
        axes[1,0].semilogy(y_axis, abs_el[el]["fluorescence"][depth,:,int((Nz-1)/2)]  , label="fluorescence" )
        axes[1,0].semilogy(y_axis, abs_el[el]["absorption"][depth,:,int((Nz-1)/2)]  , label="absorption" )
        
        axes[0,1].set_title("Dose rate J/kg/s")
        axes[0,1].plot(y_axis, abs_el[el]["dose"][depth,:,int((Nz-1)/2)]  , label="dose in the sample" )
        axes[0,1].plot(y_axis, abs_el[el]["dose_col"][depth,:,int((Nz-1)/2)]  , label="dose in the collagen" )
        axes[0,1].legend()

        #axes[1,1].set_title("Dose")
        axes[1,1].semilogy(y_axis, abs_el[el]["dose"][depth,:,int((Nz-1)/2)]  , label="dose in the sample" )
        axes[1,1].semilogy(y_axis, abs_el[el]["dose_col"][depth,:,int((Nz-1)/2)]  , label="dose in the collagen" )
        
        
        plt.savefig("{}/{}_el_{}_{:04d}.png".format(datadir,appname, el, E))
        plt.close()
        
    fig, (ax1,ax2) = plt.subplots(1,2)
    fig.suptitle("I0 absorption front E={}".format(E))
    
    ax1.set_title('Absorption')
    ax2.set_title('Dose')
    
    plt.rcParams.update({"text.usetex": True})
    
    
    
    ax1.plot(y_axis, absorption[0,:,int((Nz-1)/2)], label="Abs" ) # converting it in mu

    
    ax2.semilogy( y_axis, dose[depth,:,int((Nz-1)/2)], label="sample") # converting it in mu
    ax2.semilogy( y_axis, dose_HAP[depth,:,int((Nz-1)/2)], label="HAP") # converting it in mu
    ax2.semilogy( y_axis, dose_Col[depth,:,int((Nz-1)/2)], label="collagen from I0") # converting it in mu
    ax2.semilogy( y_axis, Tot_dose_Col[depth,:,int((Nz-1)/2)], label="collagen Total") # converting it in mu
    ax2.semilogy( y_axis, Tot_dose[depth,:,int((Nz-1)/2)], label="Total Dose") # converting it in mu
    ax2.semilogy( y_axis, abs_el["Ca"]["dose_col"][depth,:,int((Nz-1)/2)], label="Dose from Ca fluorescence") # converting it in mu
    ax2.legend()
    plt.savefig("{}/{}_Tot_{:04d}.png".format(datadir/appname, E))
    plt.close()
