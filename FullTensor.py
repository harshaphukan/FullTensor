#!/usr/bin/env python
#Class to calculate Full Strain Tensor from Energy Scan Data @ APS ID34E
#

import numpy as np
from scipy import optimize
import os


class FullTensor:
  """Initialize Deviatoric Strain Tensor (devTens), 
     Reciprocal lattice Vector (Q0), 
     Magnitude of Measured Q from Energy Scan(QMeas)
     Indices corresponding to the Q Measurement (h,k,l)
    
     Please Note: If the individual components of the R-Lattice Vectors are in row format, please convert them to Column format by taking
     the transpose.
  """

  def __init__(self,devTens,Q00,QMeas,h,k,l):

    self.devTens=devTens
    self.Q00=Q00
    self.QMeas=QMeas
    self.h,self.k,self.l=h,k,l
  
    
    
    
    
    
    
    
    
    
         
  def RealSpace(self):
    """ Given a set of reciprocal lattice vectors, this method converts them to real space vectors """    
    aStar0=np.transpose(self.Q00[:,0])
    bStar0=np.transpose(self.Q00[:,1])
    cStar0=np.transpose(self.Q00[:,2])
    


    c0=np.cross(aStar0,bStar0)
    c0=c0/np.linalg.norm(c0)
    a0=np.cross(bStar0,cStar0)
    a0=a0/np.linalg.norm(a0)
    b0=np.cross(cStar0,aStar0)
    b0=b0/np.linalg.norm(b0)
    

    
    Real00=np.vstack((a0,b0,c0))
    Real00=np.transpose(Real00)    

    return Real00  
    
    
  def DefGrad(self,Real00,Real01):
    """ Method to calculate Deformation gradient, given initial and final real space lattice vectors. """
    
    F=np.dot(Real01,np.linalg.inv(Real00))
    
    return F
     
    
  def Function(self,x):
   
    D=np.eye(3,3)+self.devTens+np.array([[x,0.,0.],[0.,x,0.],[0.,0.,x]])
    D_Inv=np.linalg.inv(D)
    
    aStar0=self.h*self.Q00[:,0]
    bStar0=self.k*self.Q00[:,1]
    cStar0=self.l*self.Q00[:,2]
    Qtheo=aStar0+bStar0+cStar0


    RHS=np.dot(D_Inv,Qtheo)
  
    Mag_RHS=np.linalg.norm(RHS)
  
  
    return np.abs(Mag_RHS-self.QMeas)
   
  def Optimize(self,Tol,UL,LL,Iter):
   """ This method uses the Brent Optimization Scheme to estimate the Hydrostatic Strain Component """
   Hydro_Opt=optimize.brent(self.Function,brack=(LL,UL),tol=Tol,full_output=True,maxiter=Iter)
   return Hydro_Opt


  



    
  
