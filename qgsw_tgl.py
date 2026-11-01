import numpy as np
from scipy.interpolate import griddata
import modgrid
import moddyn
import modelliptic

def qgsw_tgl(Htraj=None, dHi=None, c=None, lon=None, lat=None, tint=None, dtout=None, dt=None,obsspace=None,rappel=None,snu=None):
  """Tangent-linear model of qgsw() around the trajectory Htraj.

  Htraj is the output of qgsw(), i.e. one frame every dtout, while the model
  advances every dt. The linearisation point is therefore piecewise constant
  over each output window, and the accuracy of the tangent-linear model is
  limited by dtout, not by dt. Measured on ssh.nc (tint=2 days, dt=600 s) as
  the ratio of the nonlinear perturbation to the tangent-linear prediction:

      dtout      frames    |ratio - 1|
      86400 s         3      7.4e-02
      21600 s         9      9.2e-03
       3600 s        49      4.7e-04
        600 s       289      4.6e-06   (= dt, fully resolved)

  So pass dtout=dt if you need an accurate tangent-linear model, and expect
  several percent error at daily output. Note that the accuracy of qgsw_adj is
  bounded the same way.
  """
 
  way=np.sign(tint)

  ##############
  # Setups
  ##############
  dHi=dHi+Htraj[0,:,:]*0.

  grd=modgrid.grid(dHi,c,snu,lon,lat)
  time_abs=0.
  index_time=0  
  if obsspace is not None:
    dhg=np.empty((np.shape(obsspace)[0]))
    dhg[:]=np.nan
    iobs=np.where((way*obsspace[:,2]>=time_abs-dt/2) & (way*obsspace[:,2]<time_abs+dt/2))
    if np.size(iobs)>0:
      # NOTE: was `h.ravel()`, which is undefined here. At time_abs=0 the
      # perturbation field is dHi.
      dhg[iobs]=griddata((lon.ravel(), lat.ravel()), dHi.ravel(), (obsspace[iobs,0].squeeze(), obsspace[iobs,1].squeeze()))
  else:
    dhg=None

  nindex_time=int(np.abs(tint)/dtout) + 1
  dSSH=np.empty((nindex_time,grd.ny,grd.nx))
  dSSH[index_time,:,:]=dHi  

  nstep=int(abs(tint)/dt)
  stepout=int(dtout/dt)


 
  ############################
  # Active variable initializations
  ############################

  dh=+dHi
  dq,=modelliptic.h2pv(dh,grd)

  dhb=+dh # just for hguess
 
  ############################
  # Time loop
  ############################

  for step in range(nstep): 
    #print step
    time_abs=(step+1)*dt
    if (np.mod(step+1,stepout)==0):
      index_time += 1

    ############################
    #Tangent update on current trajectory
    ############################

    # Linearise around the most recent stored trajectory frame at or BEFORE the
    # current time. index_time has already been incremented above for the output
    # slot, so using it here reads one frame ahead on the last step of every
    # output window, which biased the tangent-linear model by several percent.
    # Htraj only holds one frame per dtout, so step//stepout is the best
    # piecewise-constant linearisation point available.
    h=Htraj[step//stepout,:,:].squeeze()
    q,=modelliptic.h2pv(h,grd)
    u,v, = moddyn.h2uv(h,grd)

    ############################
    # Initializations
    ############################

    dhguess=2*dh-dhb
    dhb=+dh
    dqb=+dq

    ########################
    # Main routines
    ########################

    # 1/
    du,dv, = moddyn.h2uv(dhb,grd)

    # 2/
    drq, = moddyn.qrhs_tgl(du,dv,dqb,u,v,q,grd,way)

    # 3/
    if rappel is not None:    
      dq = dqb + dt*(drq-rappel*(dqb))
    else:
      dq =dqb + dt*drq

    # 4/ From new q, we update h
    dh,=modelliptic.pv2h(dq,dhguess,grd)

    ############################
    #Saving outputs
    ############################

    if (np.mod(step+1,stepout)==0): 
      dSSH[index_time,:,:]=dh

    if obsspace is not None:
      iobs=np.where((way*obsspace[:,2]>=time_abs-dt/2) & (way*obsspace[:,2]<time_abs+dt/2))
      if np.size(iobs)>0:
        dhg[iobs]=griddata((lon.ravel(), lat.ravel()), dh.ravel(), (obsspace[iobs,0].squeeze(), obsspace[iobs,1].squeeze()))

  return dSSH,dhg




