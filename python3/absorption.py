#!/usr/bin/env python3

__author__ = "Ivo Zizak (ivo.zizak@helmholtz-berlin.de)"
__version__ = "$ 0.09 $"
__date__ = "$Date: 2021/12/23 00:00:00 $"
__copyright__ = "Copyright (c) 2021 Ivo Zizak"
__license__ = "Python"


# Version 008: using number concentrations now
import sys
import numpy as np

from scipy.ndimage import gaussian_filter
from scipy import ndimage   # for convolve
from scipy import signal   # for convolve
import pandas as pd
import re
# this is needed only for testing:
import matplotlib.pyplot as plt
import h5py
import time


from parameters import *


# just joining the file names...
appname="{}/{}".format(datadir,appname)
simappname1="{}/{}".format(datadir,simappname1)
simappname2="{}/{}".format(datadir,simappname2)


#atomic_constants=np.genfromtxt("atomic_constants.dat")
elementsList = (
    'none', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
    'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
    'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
    'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
    'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U')

#print atomic_constants[10]

        
class beam():
    def __init__():
        self.energy=0 # eV
        # intensity should be a grid
        # but this makes sense only for the parallel beam
        self.flux=0  # photon/m^2
        self.position=np.array([0.0,0.0,0.0])
        self.direction=np.array([1.0,0.0,0.0])

class material():
    def __init__(self, formula, density=None, cell_volume=None):
        self.iselement=False # this is only becauswe I can use another database to get the density
        self.ismixture=False
        
        self.formula=formula
        
        self.materials=[]
        self.elements=[]
        
        
        self.elements=self.parse_formula(formula)

        self.MolMass=self.mol_mass()
        self.Natoms=0
        print("calculated number of atoms in the cell:",self.Natoms)
        for el in self.elements:
            self.Natoms=self.Natoms+int(el[1])

        # this is just the initialisation, these depend on each other:
        self.density=density
        self.number_density=0 # how many primitive cells/molecules per m^3, unit 1/m^3
        self.set_unit_volume
        if density:
            self.set_density(density)
        if cell_volume:
            self.set_unit_volume(volume)
        
        
        print("Density=",self.number_density," [cells/m^3]")
        
        # this should be a pandas table wit diffferent columns and 
        # 1st axis Energy- 
        # I can create it when I get the energy first
        # TODO an alternative would be to load the raw tables of f1 and f2 for all elements
        # end the energy axis from the file, and interpolate later when I need the energy
        # This would make possible to change the energies without loading the file.
        # I can do this later
        self.have_table=False
                
        # absorption is atomic absorption cross section, for individual lements, 
        # as well as the total cross section.
        # should get pandas here
        #self.sigma_abs_el=self.load_absorption()
        
        # internals
        # do we have an nk-file for the material?
        # should we calculate it from the components? For now, only this is plausible, 
        # since there are not so many materials with n,k in the database.

        

    def load_absorption(self,Energy):
        # opens the elements file and loads the scattering factors
        # energy is a vector for the energies
        # keeps the individual elements absorption - for anomalous 
        #if(self.have_table==False):
            #self.table=pd.DataFrame (Energy,columns = ['Energy'] )
        self.table=pd.DataFrame (index=Energy )
        self.have_table=True
        
        r=h5py.File("Elements_01.hdf5", "r")
        print("Reading file Elements")
        self.table["mu"]=0*Energy
        lmbda=12.39841662513396/(Energy/1000.0)*1e-10 # in Angstrom to meter
        r_e = 2.8179403E-15  # in m

        self.mol_mass()

        for el in self.elements:
            print("Adding column f2_{}".format(el[0])),
            x=np.array(r["/Henke_f1f2/{}".format(el[0])])
            N=int(el[1])
            f2=np.interp(Energy, x[0], x[2])
            self.table["f2_{}".format(el[0])]=f2
            #print(self.table)
            
            # here calculate the linear absorption factor
            # now convert these into the absorption cross section
            #I0=1 # photon/m^2 - flux, beam density
    
            # From the density calculate the LINEAR absorption constant,
            # as well as the portions for every element
            # This is already multiplied by N, so it absorption of all atoms of the same species
            self.table["mu_{}".format(el[0])]=2*r_e*lmbda*( N * f2)* self.number_density
            print("Adding column mu_{}".format(el[0]))
            mass=r["/AtomicProperties/{}/AtomicMass".format(el[0])][0]
            #print("mass",type(mass))
            #print("mass",mass)
            #print("MolMass",type(self.MolMass))
            W_frac = N*mass/self.MolMass
            
            N_frac=N/self.Natoms
            
            print("Element",el[0])
            print("   mass fracion:",W_frac)
            print("    num fracion:",N_frac)
            
            self.table["mu"]=self.table["mu"]+self.table["mu_{}".format(el[0])] # /N_frac*W_frac
            
            
            #self.table["mu"]=self.table["mu"]+2*r_e*lmbda*( N*f2)* self.number_density
        if(False):
            # this is to check the calculation precision
            # it worked once, so it is not important,
            # the precision is better than 1E-16
            print("checking the precision")
            a=self.table["mu"]
            b=a
            for key in self.table.keys():
                if key.startswith("mu_"):
                    a=a-self.table[key]
            print(a/b)
        r.close()
        
        
    def parse_formula(self,formula):
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        if len(elements)==1:
            self.iselement=True # this is only becauswe I can use another database to get the density
        # run over the first part of the touple and check if this is an element.
        # if the second part is "" (no number, set it to 1, else evaluate the string to get the number
        return elements

    def mol_mass(self):
        r=h5py.File("Elements_01.hdf5", "r")
        self.MolMass=0
        for el in self.elements:
            mass=r["/AtomicProperties/{}/AtomicMass".format(el[0])]
            print(el[0],mass[0])
            self.MolMass=self.MolMass+float(el[1])*float(mass[0])
        r.close()
        print("calculated unit Mass",self.MolMass)
        return self.MolMass
    
    def set_unit_volume(self,vol):
        self.number_density=1/vol
    
    def set_number_density(self,dens):
        # number_density is in molecules (or unit cells, corresponding to the formula) per m^3
        self.number_density=dens
        AMU=1.66053907E-27# kg
        self.density=(dens/1000)*self.MolMass*AMU
        
    def set_density(self,density):
        #density is in g/cm^3
        self.density=density
        self.mol_mass()
        AMU=1.66053907E-27# kg
        # density in kg/m^3
        self.number_density=(density*1000)/self.MolMass/AMU

        
    def get_property(self,energy):
        pass


def export_as_vtk2(filename,array,origin=[0,0,0], spacing=[1,1,1]):
    f = open(filename, "w")
    f.write("# vtk DataFile Version 2.0\n")
    f.write("A Simple 3D Array of values\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_POINTS\n")
    f.write("DIMENSIONS {} {} {}\n".format(array.shape[0],array.shape[1],array.shape[2]))
    f.write("ORIGIN {} {} {}\n".format(origin[0],origin[1],origin[2]))
    f.write("SPACING {} {} {}\n".format(spacing[0],spacing[0],spacing[0]))
    f.write("POINT_DATA {}\n".format(array.shape[0]*array.shape[1]*array.shape[2]))
    f.write("SCALARS values float\n")
    f.write("LOOKUP_TABLE default\n")

    for z in range(array.shape[2]):
        for y in range(array.shape[1]):
            for x in range(array.shape[0]):
                f.write("{} ".format(array[x,y,z]))
        
    f.close()





# calculate hte fluorescence distribution here, outside of the loop
# this takes a lot of memory, so I have to reuse the arrays.
def calc_fluo_absorption( name, mu, Vsize=1e-6, size=(100,100,100) ):
    """
    First I would want to calculate only the fluorescence distribution. in a 3D array
    
    This distribution would be used only to calculate the absorption in the block. However, 
    the distribution of the fluorescence sources depends on energy, so it would have to be recalculated every time.
    Here I would only calculate the fluorescence distribution and do the convolution later.
    
    mu and Vsize are both in meter
    
    Arguments: material (I actually need only mu, size of the voxel, Nx, Ny, Nz 
    returns np.array of the shape (2*Nx+1, 2*Ny+1, 2*Nz+1) where the middle voxel is a uniform source of the fluorescence
    Values are the portions of the flourescence from the middel pixel absorbed in every single pixel. 
    If the array is large enough, teh sum should be 1
    
    TODO: it is complicated, but it keeps the arrayy smaller: it can be symmetrized.
    
    """


    # compared to ray simulations, this calculation is OK for the xoxels far away from the source, 
    # but it overestimates the "0" voxel and underestimetes the first few neighbours.
    # however, the errors are around 10%
    fluorescence_E={
        "C":277,
        "N":392,
        "O":524,
        "P":2010,
        "Ca":3690
    }
 
    Nx,Ny,Nz = size

    filename="Fluo_dist_{}".format(name)
    filename = "{}_Fluo_dist_{:04d}keV.npz".format(simappname1,fluorescence_E[name])

    try:
        saved = np.load(filename)
        FL_dist=saved["dist"]
        
        print("Loaded {} .".format(filename) )
        Sx,Sy,Sz= FL_dist.shape
        if (Sx==(2*Nx+1)) and (Sy==(2*Ny+1)) and (Sz==(2*Nz+1)):
            return FL_dist
            
        print("found the distribution, but sizes mismatch")
        # recalculate
        filename="{}_Fluo_dist_{}.npz".format(simappname2,name)
    except:
        print("could not use {} , trying the other one".format(filename) )

        filename="{}_Fluo_dist_{}.npz".format(simappname2,name)
        print("Unexpected error:", sys.exc_info()[0])
        # if something went wrong, just recalculate it
        try:
            saved = np.load(filename)
            FL_dist=saved["dist"]
            
            print("Loaded {} .".format(filename) )
            Sx,Sy,Sz= FL_dist.shape
            if (Sx==(2*Nx+1)) and (Sy==(2*Ny+1)) and (Sz==(2*Nz+1)):
                return FL_dist
        except:
            print("Unexpected error:", sys.exc_info()[0])
            

    print("could not use {} , recalculating".format(filename) )
    Fluo_dist=np.ndarray([2*Nx+1,2*Ny+1,2*Nz+1])
    FL_dist=np.ndarray([2*Nx+1,2*Ny+1,2*Nz+1])
    for dx1 in np.arange(2*Nx+1):
        for dy1 in np.arange(2*Ny+1):
            for dz1 in np.arange(2*Nz+1):
                r12s=((dx1-Nx)*(dx1-Nx) + (dy1-Ny)*(dy1-Ny) + (dz1-Nz)*(dz1-Nz))*Vsize*Vsize
                r12=np.sqrt(r12s)
                if 2*r12s>Vsize*Vsize:
                    SolidAngle= Vsize * Vsize /(4* np.pi * r12s) 
                    FL_dist[dx1,dy1,dz1] = SolidAngle * (np.exp(-1*mu*r12)-np.exp(-1*mu*(r12+Vsize)) )
                else:
                    # the same voxel absorption
                    FL_dist[dx1,dy1,dz1] = (1-np.exp(-1*mu*(Vsize/2)))

    np.savez_compressed(filename, dist=FL_dist, mu=mu, Vsize=Vsize, )


    # normalize the sum to 1 (experimental)
    sumFL = FL_dist.sum()
    FL_dist = FL_dist / sumFL

    return FL_dist

# I hold all the data for elements in the database
# Elements_01.hdf5
# when I need some data for this element, I can call the function to get it grom file
# It does not make sense to load all the data at the begin of the program and keep it there to the end, there is really a lot of data in the database.

# test database:
def test():
    # read test
    print("Generating HAP")
    HAP = material("Ca10P6O26H2",density=3.18) #density in g/cm^3
    Col = material("C65H102N18O21",density=1.35)
    
  

    E_f=[]
    for el in fluorescence_E.keys():
        E_f.append(fluorescence_E[el])
    #E_f=[524, 2010, 2135, 3690, 4013]
    E_f.extend(E_beam)
    Energy=np.array(E_f)
    #Energy=np.array[3600,8000,10000,18000]

    print("All energies:", Energy, "for the tables")

    # The primary beam has the continuous energy distribution, 
    # so the primary beam intensity is given by I0
    # I0 is distributed across the beam cross section,
    # The flux is Phi0 = I0 / beam_cross_section_area

    # filling the absorption columns in the maaterial properties with the
    # data from the hp5 filfe
    HAP.load_absorption(Energy)
    Col.load_absorption(Energy)
    
    # now convert these into the absorption cross section
    # calculating the linear absorption of HAP
    # lmbda=12.39841662513396/(Energy/1000.0)*1e-10 # in Angstrom to meter
    r_e = 2.8179403E-15  # in m
    
    
    #########################################################################################################################
    # let us calculate the fluorescence from the Ca
    # for this we need a model, we can't just calculate the intensities like this
    
    # Model: 1 mu blocks in the beam, every consisting of a mixture of HAP and Collagen
    # 10 mu beam size, 1 mm sample thickness
    # Probability that the initial photon is absorbed by a calcium atom in mu_HAP
    # we calculate the transmission of 1 mu Ca in 100% HAP sample

    # TODO here I removed the parameter definitions!
    
    #let us make a beam profile
    Nx=int(sample_thickness/Vsize)
    Ny=int(2*box_Y_size/Vsize)
    Nz=int(2*box_Z_size/Vsize)
    ycenter=int(Ny/2)
    zcenter=int(Nz/2)

    x_axis=np.linspace(0, sample_thickness, num=Nx)
    y_axis=np.linspace(-box_Y_size, box_Y_size, num=Ny)
    z_axis=np.linspace(-box_Z_size, box_Z_size, num=Nz)

    zs,ys=np.meshgrid(z_axis, y_axis) # ys, like on surface
    
    # circle
    beam_profile=np.where(ys*ys+zs*zs<beam_dia*beam_dia/4.0,1.0, 0.0)
    # gaussian blur
    #TODO release this at the end
    GaussBlurSigma=2
    beam_profile=ndimage.gaussian_filter(beam_profile, sigma=GaussBlurSigma)
    
    # top hat gauss
    #beam_profile=(np.exp(-1*np.power((ys*ys+zs*zs)/(beam_dia*beam_dia),5)))
    
    # noemalize
    IN=beam_profile.sum()
    beam_profile=beam_profile/IN*I0 # norming the beam profile to the total intensity of I0

    print("shape of y_axis:",y_axis.shape)
    print("shape of z_axis:",z_axis.shape)
    print("shape of ys:",ys.shape)
    print("shape of zs:",zs.shape)
    print("shape of beam_profile:",beam_profile.shape)
    
    print(beam_profile[zcenter,:])

    # calculate the FWHM
    knife_scan=np.sum(beam_profile,axis=1)
    axis=y_axis
    
    iMAX=np.argmax(knife_scan)
    xMAX=axis[iMAX]
    axis=axis-xMAX
    MAX=knife_scan[iMAX]
    

    C=np.where(knife_scan>MAX/2)
    # indices of the first and last value larger than MAX/2
    left=C[0][0]
    right=C[0][-1]
    #print(C[0])
    
    # interpolate left and right to find the x-value where the line crosses MAX/2
    xleft= axis[left-1]+(MAX/2-knife_scan[left-1]) * (axis[left]-axis[left-1])/(knife_scan[left]-knife_scan[left-1])
    xright=axis[right]+(MAX/2-knife_scan[right]) * (axis[right+1]-axis[right])/(knife_scan[right+1]-knife_scan[right])
    print("FWHM={:4f}".format((xright-xleft)*1E6))


    knife_scan=np.sum(beam_profile,axis=1)

    # beam_profile sets how many photons arrive into one voxel at the surface

    # our sample is a mixture of two materials, conc_HAP is VOLUME portion of HAP
    # ALL voxels in this calculation have the same composition (on micrometer range)
    conc_HAP=0.4
    conc_Col=1.0-conc_HAP

    # in a mixture, the absorption is a sum of absorptions of the components:
    #  abs=abs_Col*conc_Col+abs_HAP*conc_HAP
    # since actually the density of HAP and Col change when they are mixed, we could have started with 
    # dinsities multiplied wit the composition brom the beginning...
    
    # So if a voxel is alone in the world, the number of photons which are goinng to be absorbed in it is:
    # I_abs=I0(1-exp(-1*mu*Vsize))
    # for every voxels hiding behind N voxels:
    # I_abs_N=I0*(exp(-N*mu*bsize) - exp(-(N+1)*mu*bsize ) ) )

    for E in E_beam:
        #Construct the absorption depth profile array
        #Weight fraction:
        W = conc_HAP * HAP.density + conc_Col * Col.density
        W_HAP = (conc_HAP * HAP.density)/W
        W_Col = (conc_Col * Col.density)/W

        # Number fraction
        N_HAP = conc_HAP * HAP.number_density
        N_Col = conc_Col * Col.number_density
        N_tot = N_HAP + N_Col
        N_HAP = N_HAP / N_tot
        N_Col = N_Col / N_tot

        print("Weight fractions")
        print("HAP:",W_HAP)
        print("Col:",W_Col)
        print("sum:",W_Col+W_HAP)

        print("Number fractions")
        print("HAP:",N_HAP)
        print("Col:",N_Col)
        print("sum:",N_Col+N_HAP)

        mu=conc_HAP*HAP.table["mu"]+conc_Col*Col.table["mu"]
        abs_fac=mu
        inc_I = np.exp(-1*mu[E]*x_axis)
        trans_I = np.exp(-1*mu[E]*(x_axis+Vsize))
        abs_I=inc_I-trans_I
        print(trans_I[0:5])
        
        # calculating how much is absorbed by specific components
        # I only need the absorption in one voxel, because I know how much photons are entering this voxel
        
        mu= conc_Col * Col.table["mu"]
        abs_fac_Col = mu

        abs_I_Col = abs_I * conc_Col * abs_fac_Col[E] / abs_fac[E]

        mu= conc_HAP * HAP.table["mu"]
        abs_fac_HAP = mu
        #print("wrong HAP:", inc_I*(1-np.exp(-1*mu[E]*(Vsize))))
        abs_I_HAP = abs_I * conc_HAP * abs_fac_HAP[E] / abs_fac[E]

        abs_I0_elem={}
        abs_I_total=0

        for elem in elements:
            print("Calculating {} absorption of I0".format(elem))
            tableT="mu_{}".format(elem)
            try:
                mu=conc_Col * Col.table[tableT]
                #abs_I0_elem_Col=inc_I*(1-np.exp(-1*mu[E]*(Vsize)))
                abs_I0_elem_Col=abs_I*conc_Col * mu[E] / abs_fac[E]
            except:
                print(" no {} in Collagen".format(elem))
                abs_I0_elem_Col=inc_I*0

            try:
                mu=conc_HAP * HAP.table[tableT]
                abs_I0_elem_HAP=abs_I * conc_HAP * mu[E] / abs_fac[E]
            except:
                print(" no {} in HAP".format(elem))
                abs_I0_elem_HAP=inc_I*0

            abs_I0_elem[elem]=abs_I0_elem_HAP+abs_I0_elem_Col
            abs_I_total=abs_I_total + abs_I0_elem_HAP+abs_I0_elem_Col

            print("percentage absorbed by", elem,":",abs_I0_elem[elem])
    
        
        print("           Total absorbed",abs_I_total )

        starttime=time.time()

        # I know we have a cube, and that the photons are comming under funny angle, but in
        # the first approximation, I only whant to concider a volume and 
        # cross section, so I assume the cube is oriented with one side to the fluorescence source

        # F3D_I0_abs is the absorption of the primary beam in one voxel
        F3D_I0_abs=np.ndarray([Nx,Ny,Nz])
        for ix in range(F3D_I0_abs.shape[0]):
            F3D_I0_abs[ix,:,:]=beam_profile*abs_I[ix]

        # this is the part of the dose from the primary beam
        #weight is weight of a voxel
        # conc_ is volume portion, and W_ is the weight portion
        weight=(conc_Col*1.35 + conc_HAP * 3.18)*Vsize*Vsize*Vsize * 1000 # g/cm^3 in kg/m^3
        dose  = F3D_I0_abs      *  E * 1.602176634E-19  / weight # convert in Joule!
        # for a given energy, the portion of the intensity absorbed by Colagen is always the same

        # fract_Col is the fraction of all photons absorbed by the collagen
        #fract_Col = 1 / ( 1 + (1-np.exp(-1*abs_fac_Col[E])) / (1-np.exp(-1*abs_fac_HAP[E]) ) )
        fract_Col = conc_Col*Col.table["mu"][E] / abs_fac[E]
        weight_Col = conc_Col * 1.35 * Vsize*Vsize*Vsize * 1000
        dose_Col  = F3D_I0_abs * fract_Col * E * 1.602176634E-19 / weight_Col

        weight_HAP = conc_HAP * 3.18 * Vsize*Vsize*Vsize * 1000
        
        #fract_HAP = 1 / ( 1 + (1-np.exp(-1*abs_fac_HAP[E])) / (1-np.exp(-1*abs_fac_Col[E]) ) )
        fract_HAP = conc_HAP*HAP.table["mu"][E] / abs_fac[E]
        
        dose_HAP  = F3D_I0_abs * fract_HAP * E * 1.602176634E-19 / weight_HAP

        print("For Energy {}".format(E))
        print("   fraction absorbed by Collagen = {}".format(fract_Col))
        print("   fraction absorbed by HAP = {}".format(fract_HAP))
        

        filename="{}_F3D_I0_{:04d}".format(appname,E)
        np.savez_compressed(filename, absorption=F3D_I0_abs, dose=dose, dose_Col=dose_Col, dose_HAP=dose_HAP )


        for elem in elements:
            # calculate the absorption of the element
            # and multiply it at the same time with the fluorescence yield to obtain the fluorescence,
            # and put it int the 3D model
            E_f=fluorescence_E[elem]

            F3D_Fluo_elem = F3D_I0_abs*0
            for ix in range(F3D_Fluo_elem.shape[0]):
                F3D_Fluo_elem[ix,:,:] = beam_profile * abs_I0_elem[elem][ix]*fluorescence_yield_K[elem]
            
            # calculate how the flourescence is absorbed
            # WARNING I have to use the complete abs_fac here, later I can calculate the Col part.
            mu=abs_fac[E_f]
            newshape=(150,150,150)

            print("Calculating {} fluorescence".format(elem))
            Fluo_dist = calc_fluo_absorption(elem,mu,Vsize,newshape)
            print("  max",Fluo_dist.max())
            print("  sum",Fluo_dist.sum())

            # F3D_Fluo_Abs is the absorbed number of photons from the fluorescence 
            # of the specified element in one voxel.
            # BUG: I should calculate this only once, like convolution of the I0 with the fluorescence, and than only multiply with the absorption...
            F3D_Fluo_Abs = signal.fftconvolve( F3D_Fluo_elem , Fluo_dist, mode='same')
            
            #fract_Col = 1 / ( 1 + (1-np.exp(-1*abs_fac_Col[E_f])) / (1-np.exp(-1*abs_fac_HAP[E_f]) ) )
            fract_Col = conc_Col*Col.table["mu"][E_f] / abs_fac[E_f]



            print("totals of absorption should be close to the total fluorescence:")
            # this is not neccessray, but it is for me to check that the convolution went OK
            print(elem,":")
            print("          total fluorescence {:g} photons".format(np.sum(F3D_Fluo_elem)))
            print(" Absorbed total fluorescence {:g} photons".format(np.sum(F3D_Fluo_Abs)))

            # DOSE: absorbed energy/mass
            # We know the density, so we can calculate the voxel weight
            
            # part absorbed in Colagen
            dose_fluo_Col = F3D_Fluo_Abs * fract_Col * E_f * 1.602176634E-19/weight_Col # convert in Joule!
            
            dose_fluo = F3D_Fluo_Abs * E_f * 1.602176634E-19/weight # convert in Joule!
            
            # for a given energy, the portion of the intensity absorbed by Colagen is always the same
            
            dose     = dose     + dose_fluo
            dose_Col = dose_Col + dose_fluo_Col
            
            # save all the results for the later plotting and analysis:
            filename="{}_F3D_{}_{:04d}".format(appname,elem,E)
            np.savez_compressed(filename, absorption = F3D_Fluo_Abs, fluorescence = F3D_Fluo_elem, dose=dose_fluo, dose_col=dose_fluo_Col )

        filename="{}_F3D_Tot_{:04d}".format(appname,E)
        np.savez_compressed(filename, absorption=F3D_I0_abs, dose=dose, dose_Col=dose_Col )


test()
