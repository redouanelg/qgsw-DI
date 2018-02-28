import numpy
from math import cos,sin,pi,isnan
import matplotlib.pylab as plt
import pdb

def h2pv(h,grd):
    """ SSH to Q

    Args:
        h (2D array): SSH field.
        grd (Grid() object): check modgrid.py

    Returns:
        q: Potential Vorticity field  
    """
    g=grd.g
    ny,nx,=numpy.shape(grd.mask)
    q=- g*grd.f0/(grd.c**2) *h

    q[1:-1,1:-1]=g/grd.f0[1:-1,1:-1]*((h[2:,1:-1]+h[:-2,1:-1]-2*h[1:-1,1:-1])/grd.dy[1:-1,1:-1]**2 \
                                      + (h[1:-1,2:]+h[1:-1,:-2]-2*h[1:-1,1:-1])/grd.dx[1:-1,1:-1]**2) \
    - g*grd.f0[1:-1,1:-1]/(grd.c[1:-1,1:-1]**2) *h[1:-1,1:-1]
    # +grd.f0[1:-1,1:-1] - g*grd.f0[1:-1,1:-1]/(grd.c[1:-1,1:-1]**2) *h[1:-1,1:-1]
    ind=numpy.where((grd.mask==1))
    q[ind]=- g*grd.f0[ind]/(grd.c[ind]**2) *h[ind]
    #qq[ind]=grd.f0[ind] - g*grd.f0[ind]/(grd.c[ind]**2) *h[ind]
    #qtemp=- g*grd.f0/(grd.c**2) *h
    #q[numpy.where((numpy.isnan(q)))]=qtemp[numpy.where((numpy.isnan(q)))]
    ind=numpy.where((grd.mask==0))
    q[ind]=0

    return q,


def get_aq_from_ah(ah,grd):

    g=grd.g
    ny,nx,=numpy.shape(grd.mask)
    aq=numpy.empty((ny,nx,2))
    qq= - g*grd.f0/(grd.c**2) *h

    qq[1:-1,1:-1]=g/grd.f0[1:-1,1:-1]*((h[2:,1:-1]+h[:-2,1:-1]-2*h[1:-1,1:-1])/grd.dy[1:-1,1:-1]**2 \
                                       + (h[1:-1,2:]+h[1:-1,:-2]-2*h[1:-1,1:-1])/grd.dx[1:-1,1:-1]**2) \
    - g*grd.f0[1:-1,1:-1]/(grd.c[1:-1,1:-1]**2) *h[1:-1,1:-1]
    # +grd.f0[1:-1,1:-1] - g*grd.f0[1:-1,1:-1]/(grd.c[1:-1,1:-1]**2) *h[1:-1,1:-1]
    ind=numpy.where((grd.mask==1))
    qq[ind]=- g*grd.f0[ind]/(grd.c[ind]**2) *h[ind]
    #qq[ind]=grd.f0[ind] - g*grd.f0[ind]/(grd.c[ind]**2) *h[ind]
    ind=numpy.where((grd.mask==0))
    qq[ind]=0
    q[:,:,0]=qq
    q[:,:,1]=qq

    return aq,


def pv2h(q,hg,grd,nitr=1):
    """ Q to SSH
    
    This code solve a linear system of equations using Conjugate Gradient method

    Args:
        q (2D array): Potential Vorticity field
        hg (2D array): SSH guess
        grd (Grid() object): check modgrid.py

    Returns:
        h (2D array): SSH field. 
    """
    ny,nx,=numpy.shape(hg)
    g=grd.g


    x=hg[grd.indi,grd.indj]
    q1d=q[grd.indi,grd.indj]

    #aaa=g/grd.f01d
    #bbb=-g*grd.f01d/grd.c1d**2
    #ccc=q1d-grd.f01d
    aaa=g/grd.f01d
    bbb=-g*grd.f01d/grd.c1d**2
    ccc=+q1d
    #ccc=+q1d-grd.f01d  

    aaa[grd.vp1]=0
    bbb[grd.vp1]=1
    ccc[grd.vp1]=x[grd.vp1]  ##boundary condition

    vec=+x

    avec,=compute_avec(vec,aaa,bbb,grd)
    gg=avec-ccc
    p=-gg
    #print 'test1', numpy.var(gg)

    for itr in range(nitr-1): 
        vec=+p
        avec,=compute_avec(vec,aaa,bbb,grd)
        tmp=numpy.dot(p,avec)
        
        if tmp!=0. : s=-numpy.dot(p,gg)/tmp
        else: s=1.
        
        a1=numpy.dot(gg,gg)
        x=x+s*p
        vec=+x
        avec,=compute_avec(vec,aaa,bbb,grd)
        gg=avec-ccc
        #print 'test', numpy.var(gg)
        a2=numpy.dot(gg,gg)
        
        if a1!=0: beta=a2/a1
        else: beta=1.
        
        p=-gg+beta*p
    
    vec=+p
    avec,=compute_avec(vec,aaa,bbb,grd)
    val1=-numpy.dot(p,gg)
    val2=numpy.dot(p,avec)
    if (val2==0.): 
        s=1.
    else: 
        s=val1/val2
    
    #pdb.set_trace()
    a1=numpy.dot(gg,gg)
    x=x+s*p

    # back to 2D
    h=numpy.empty((ny,nx))
    h[:,:]=numpy.NAN
    h[grd.indi,grd.indj]=x[:]


    return h,



def compute_avec(vec,aaa,bbb,grd):
    
    avec=numpy.empty(grd.np,) 
    #avec[grd.vp2]=aaa[grd.vp2]*((vec[grd.vp2e]+vec[grd.vp2w]-2*vec[grd.vp2])/(grd.dx1d[grd.vp2]**2)+(vec[grd.vp2n]+vec[grd.vp2s]-2*vec[grd.vp2])/(grd.dy1d[grd.vp2]**2)) + bbb[grd.vp2]*vec[grd.vp2]
    avec[grd.vp2]=aaa[grd.vp2]*((vec[grd.vp2e]+vec[grd.vp2w]-2*vec[grd.vp2])/(grd.dx1d[grd.vp2]**2)+(vec[grd.vp2n]+vec[grd.vp2s]-2*vec[grd.vp2])/(grd.dy1d[grd.vp2]**2)) + bbb[grd.vp2]*vec[grd.vp2]
    avec[grd.vp1]=vec[grd.vp1]
 
    return avec,


