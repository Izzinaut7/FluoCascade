import numpy as np

#BUG does not work yet!

from parameters import *

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

    csv_filename="{}/{}_F3D_{}_{:04d}.csv".format(datadir,appname,"I0",E)
    np.savetxt(csv_filename, [dose,absorption,dose_HAP,dose_Col,Tot_dose,Tot_dose_Col] )
    

# BUG Ivo @Katrein : this does not do anything, should I remove it?
#filename="{}/{}_F3D_{}_{:04d}.npz".format(datadir,appname,"I0",E)

#data = np.load()

#for key, value in data.items():
#    np.savetxt(r'/data-zmk/basic-res-shared/tomo/Katrein/3rd_paper_DAMAGE/IVO_SIMULATION/ALLE_DATEN/secondary_emission' + fluorescence + ".csv", value)
    

