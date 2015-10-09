#!/bin/python

from scipy import constants
import scipy as sp
import sys
import math as mt

def input():
 try:
  inputfile=sys.argv[1]
  ifile=open(inputfile,'r')
  while 1:
   line=ifile.readline()
   if not line: break
   if not line.startswith('#'):
    line=line.rstrip()
    exec(line)  
  ifile.close()
 except:
  print 'Usage:',sys.argv[0],'inputfile','outputfile'; sys.exit(1) 
 return Npart,float(Lbox),float(mass),float(eps),float(sigma),float(T),float(dTiSt),NTiSt,float(Rcut)

def iout(Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt,Rcut):
 hellomsg=""" MS-LJ v.1.0\n (c) 2014, Sviataslau V. Kohut
 ------------------------------------------
 INPUT DATA\n
 Parameter of the system: 
 A number of particles = %5d
 Dimension of the box = %6.3f (m*10**10)
 Mass of a particle = %6.3f (a.u.)
 Temperature = %5.2f K\n
 Integration (velocity-Verlet):
 Time step = %6.3f fs
 A number of time steps = %5d\n
 Lennard-Jones (LJ) potential:
 Epsilon = %4.1f K (in units of k_B)
 Sigma = %6.3f (m*10**10)
 Cutoff radius for interparticle interactions = %4.2f (in units of the LJ sigma)
 --------------------------------------------""" % (Npart,Lbox,mass,T,dTiSt,NTiSt,eps,sigma,Rcut)
 print hellomsg
 return

def convr(Lbox,T,dTiSt,mass,eps):
 Lbox=Lbox/sigma # length to MD
 T=T/eps # temperature to MD
 au=1.660538921e-27 # in SI
 mass=mass*au # mass of a particle in kg
 eps=eps*sp.constants.k # energy to SI
 tMD=sigma*10**(-10)*mt.sqrt(mass/eps) #MD unit of time in SI
 dTiSt=dTiSt*10**(-15)/tMD # time to MD
 return

Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt,Rcut=input() 

iout(Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt,Rcut)

convr(Lbox,T,dTiSt,mass,eps)

