import numpy as np
from math import cos,sin,pi,isnan
from scipy.interpolate import griddata
import time
import numpy.matlib as matlib
import modgrid
import moddyn
import modelliptic
import matplotlib.pylab as plt
import pdb

def qgsw(Hi=None, c=None, lon=None, lat=None, tint=None, dtout=None, dt=None,obsspace=None,Hm=None,rappel=None,snu=None):
    """ QG Shallow Water model

    Args:
        Hi (2D array): Initial SSH field.
        c (same size as Hi): Rossby first baroclinic phase speed
        lon (2D array): longitudes
        lat (2D array): latitudes
        tint (scalar): Time of propagator integration in seconds. Can be positive (forward integration) or negative (backward integration)
        dtout (scalar): Time period of outputs
        dt (scalar): Propagator time step

    Returns:
        SSH: 3D array with dimensions (timesteps, height, width), SSH forecast  
    """
    way=np.sign(tint)

  ##############
  # Setups
  ##############

    grd=modgrid.grid(Hi,c,snu,lon,lat)
    #plt.figure()
    #plt.pcolor(grd.mask)
    #plt.show()
    time_abs=0.
    index_time=0  
    if obsspace is not None:
        hg=np.empty((np.shape(obsspace)[0]))
        hg[:]=np.NAN
        iobs=np.where((way*obsspace[:,2]>=time_abs-dt/2) & (way*obsspace[:,2]<time_abs+dt/2))
        if np.size(iobs)>0:
            hg[iobs]=griddata((lon.ravel(), lat.ravel()), h.ravel(), (obsspace[iobs,0].squeeze(), obsspace[iobs,1].squeeze()))
    else:
        hg=None

    nindex_time=np.abs(tint)/dtout + 1
    SSH=np.empty((nindex_time,grd.ny,grd.nx))
    SSH[index_time,:,:]=Hi  

    nstep=int(abs(tint)/dt)
    stepout=int(dtout/dt)

  ############################
  # Passive variable initializations
  ############################
 
    if rappel is not None:
        Qm,=modelliptic.h2pv(Hm,grd)
        #Qm[np.where(grd.mask==1)]=q[np.where(grd.mask==1)]

  ############################
  # Active variable initializations
  ############################

    h=+Hi
    q,=modelliptic.h2pv(h,grd)

    hb=+h # just for hguess

  ############################
  # Time loop
  ############################

    for step in range(nstep): 
        #print step
        time_abs=(step+1)*dt
        if (np.mod(step+1,stepout)==0):
            index_time += 1

        ############################
        #Initialization of previous fields
        ############################

        hguess=2*h-hb
        hb=+h
        qb=+q    

        ########################
        # Main routines
        ########################

        # 1/ 
        u,v, = moddyn.h2uv(h,grd)

        # 2/
        rq, = moddyn.qrhs(u,v,qb,grd,way)

        # 3/    
        if rappel is not None:
            q = qb + dt*(rq-rappel*(qb-Qm))
        else:
            q =qb + dt*rq

        # 4/
        h,=modelliptic.pv2h(q,hguess,grd)


        ############################
        #Saving outputs
        ############################

        if (np.mod(step+1,stepout)==0): 
            SSH[index_time,:,:]=h

        if obsspace is not None:
            iobs=np.where((way*obsspace[:,2]>=time_abs-dt/2) & (way*obsspace[:,2]<time_abs+dt/2))
            if np.size(iobs)>0:
                hg[iobs]=griddata((lon.ravel(), lat.ravel()), h.ravel(), (obsspace[iobs,0].squeeze(), obsspace[iobs,1].squeeze()))


    return SSH,hg




