#!/usr/bin/python

#packages
import math as mt
import numpy as np
import scipy as sp
from scipy import constants
from scipy import spatial
import sys
#how to define: sp/np.X.function 

# pick up input parameters
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
# set default value for TiStPos if there is no input provided
  try:
   TiStPos
  except NameError:
   TiStPos=[]
  for i in TiStPos:
   if i < 0 or i > NTiSt:
    print 'Check your input parameters for TiStPos!'
    sys.sleep(3)
 except:
  print 'Usage:',sys.argv[0],'inputfile','>','outputfile'; sys.exit(1) 
 return Npart,float(Lbox),float(mass),float(eps),float(sigma),float(T),float(dTiSt),NTiSt,TiStPos,float(Rcut),float(Rl)

# convertion of input data:
# length, temperature,time to MD units
# energy to SI
def convr(Lbox,T,dTiSt,mass,eps):
 Lbox=Lbox/sigma # length to MD
 T=T/eps # temperature to MD
 au=1.660538921e-27 # in SI
 mass=mass*au # mass of a particle in kg
 eps=eps*sp.constants.k # energy to SI
 tMD=sigma*10**(-10)*mt.sqrt(mass/eps) #MD unit of time in SI
 dTiSt=dTiSt*10**(-15)/tMD # time to MD
 return Lbox,T,dTiSt

# summarize and return input data
def iout(Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt):
 hellomsg="""
 ***************************************************************************** 
                                    MODY-LJ v.1.0 
                          (c) 2014, Sviataslau V. Kohut
 *****************************************************************************
 INPUT DATA\n
 Parameters of the system: 
 A number of particles = %5d
 Dimension of the box = %4.2f (m*10**10) = %4.2f (m.d.u)
 Mass of a particle = %6.3f (a.u.)
 Temperature = %5.2f K\n
 Integration (velocity-Verlet):
 Time step = %6.3f fs
 A number of time steps = %5d\n
 Lennard-Jones (LJ) potential:
 Epsilon = %4.1f K (in units of k_B)
 Sigma = %6.3f (m*10**10)
 Cutoff radius for interparticle interactions = %4.2f (in units of the LJ sigma)
 Switching parameter for the potential = %4.2f (in units of the LJ sigma)
 -------------------------------------------------------------------------------""" % (Npart,Lbox,Lbox/sigma,mass,T,dTiSt,NTiSt,eps,sigma,Rcut,Rl)
 print hellomsg
 return

# generate the uniform cubic grid
def GenGrd(Npart,Lbox):
 Np=int(mt.ceil(mt.pow(Npart,1./3.))+2) # a number of points to generate in each dimension
 xGrd=np.linspace(0,Lbox,Np)
 yGrd=np.linspace(0,Lbox,Np)
 zGrd=np.linspace(0,Lbox,Np)
 Grid=np.array([],dtype=float).reshape(0,3)
 counter=0
 for k in xrange(1,Np-1):
  for j in xrange(1,Np-1): 
   for i in xrange (1,Np-1):
    counter=counter+1
    if counter > Npart:
     break
    else:            
     Grid=np.vstack([Grid,[xGrd[i],yGrd[j],zGrd[k]]])
 return Grid 

# randomly generate initial velocities for a system of particles
def GenVel(T):
 VelMod=np.random.random_sample(Npart,)
 VelMod=VelMod/np.sum(VelMod) # normalize
 Ek=1.5*Npart*T # in MD units Ekin(eps) = 3/2 * T(eps/k_B) *N
 for i in xrange(Npart):
  VelMod[i]=mt.sqrt(2*VelMod[i]*Ek) # modulus of velocity per particle 
  Vel=np.random.random_sample(3*Npart,)
 for i in xrange(3*Npart):
  Vel[i]=2*Vel[i]-1 # now it is randomly distributed in (-1,1)
 Vel=Vel.reshape(Npart,3)
 for i in xrange(Npart):
  sumsq=mt.sqrt(Vel[i][0]**2+Vel[i][1]**2 +Vel[i][2]**2)
  Vel[i]=Vel[i]/sumsq
  Vel[i]=Vel[i]*VelMod[i] # normalize
 cormom(Vel,Npart) # correction for a non-zero total momentum
 scf= mt.sqrt(Ek/Ekin(Vel)) #scaling factor
 Vel=scf*Vel 
 return Vel

# calculate the potential energy   
def Epot(Pos):
 Utot=0.
 sfijxrij=0.
 rijs=np.array([],dtype=float)
 for i in xrange(Npart):
  j=0
  while j<i:
   rij=sp.spatial.distance.euclidean(Pos[j],Pos[i])
   rijs=np.append(rijs,rij) 
   if rij > Rcut:
    j=j+1
    pass
   else:
    uij=sLJ(rij,Rcut,Rl)*uLJ(rij)
    rijv=Pos[i]-Pos[j]
    fijv=sLJ(rij,Rcut,Rl)*fLJ(rij,rijv)
    df=np.dot(fijv,rijv)
    sfijxrij=sfijxrij+df
    Utot=Utot+uij
    j=j+1
  rijs=np.sort(rijs,kind='mergesort') 
 return Utot,sfijxrij,rijs # potential energy + sum_ij fij.rij [dot product needed to calculate pressure]

# calculate the kinetic energy
def Ekin(VelXYZ):
 Energy=0.
 for i in xrange(Npart):
  velmod=mt.sqrt(VelXYZ[i][0]**2+VelXYZ[i][1]**2 +VelXYZ[i][2]**2) # in MD units!
  Energy=Energy+0.5*velmod**2
 return Energy

# evaluate the LJ potential for for ij particle interation
def uLJ(rij):
 uij = 4.*(rij**(-12)-rij**(-6)) 
 return uij

# scaling factor for the modified LJ potential
def sLJ(rij,Rcut,Rl):
 if rij <= Rl:
  S=1.
 elif rij > Rl and rij < Rcut:
  S=1.-((rij-Rl)**2*(3*Rcut-Rl-2*rij))/((Rcut-Rl)**3)
 elif rij >= Rcut:
  S=0.
 return S 

# adjust the coordinates for PBC
def PBC(Pos):
 for i in xrange(Npart):
# x-component
  if Pos[i][0] < 0.:
   Pos[i][0] = Pos[i][0] + Lbox
  if Pos[i][0] > Lbox:
   Pos[i][0] = Pos[i][0] - Lbox
# y-component
  if Pos[i][1] < 0.:
   Pos[i][1] = Pos[i][1] + Lbox
  if Pos[i][1] > Lbox:
   Pos[i][1] = Pos[i][1] - Lbox
# z-component
  if Pos[i][2] < 0.:
   Pos[i][2] = Pos[i][2] + Lbox
  if Pos[i][2] > Lbox:
   Pos[i][2] = Pos[i][2] - Lbox

# for later use... 
# calculate the force acting on ith particle
#def Force(N,Grid):
# Fi=np.zeros(3)
# for i in xrange(Npart):
#  if i!=N-1: # self-interaction
#   rij=sp.spatial.distance.euclidean(Grid[N-1],Grid[i])
#   rijv=Grid[N-1]-Grid[i]
#   if rij < Rcut:   
#    fijv=sLJ(rij,Rcut,Rl)*fLJ(rij,rijv)
#   else:   
#   fijv=fLJ(rij,rijv)
#   Fi=Fi+fijv 
# return Fi

# generate the array of forces
def Forces(Pos):
 allforces=np.zeros(3*Npart).reshape(Npart,3)
 for i in xrange(Npart):
  j=0
  while j<i:
    rij=sp.spatial.distance.euclidean(Pos[i],Pos[j])       
    rijv=Pos[i]-Pos[j]
    fijv=fLJ(rij,rijv)
    allforces[i]=allforces[i]+fijv
    allforces[j]=allforces[j]-fijv
    j=j+1
 return allforces

# calculate the interaction force 
def fLJ(rij,rijv):
 pref=48.*(rij**(-14)-0.5*rij**(-8))
 fijv=rijv*pref
 return fijv

# correct the velocities for the non-zero total momentum
def cormom(VelXYZ,Npart):
 tmom=np.sum(VelXYZ,axis=0) # returns the vector (sum px,sum py,sum pz)
 tmomx=tmom[0]
 tmomy=tmom[1]
 tmomz=tmom[2]
 if tmomx > 1.e-12 or tmomy > 1.e-12 or tmomz > 1.e-12:
  momx=tmomx/Npart
  momy=tmomy/Npart
  momz=tmomz/Npart
  for j in xrange(Npart):
   VelXYZ[j][0]=VelXYZ[j][0]-momx # x-component
   VelXYZ[j][1]=VelXYZ[j][1]-momy # y-component
   VelXYZ[j][2]=VelXYZ[j][2]-momz # z-component

# Velocity-Verlet algorithm
def VelVer(Pos,InitPos,VelXYZ,dTiSt,NTiSt):
# save some data to analyze later
 KinEgs=np.empty(NTiSt) # kinetic energies
 allFixRi=np.empty(NTiSt) # sum of F_ . R_i (pressure)
 allvar=np.empty(NTiSt)  # variance (diffusion)
 for i in xrange(NTiSt):
# calculate forces for initial positions
  fatt=Forces(Pos)
# copy initial positions  
  InitPos=np.copy(Pos)
  for j in xrange(Npart):
# update positions   
   Pos[j][0]=Pos[j][0]+VelXYZ[j][0]*dTiSt+0.5*fatt[j][0]*dTiSt**2 # x-component
   Pos[j][1]=Pos[j][1]+VelXYZ[j][1]*dTiSt+0.5*fatt[j][1]*dTiSt**2 # y-component 
   Pos[j][2]=Pos[j][2]+VelXYZ[j][2]*dTiSt+0.5*fatt[j][2]*dTiSt**2 # z-component
# periodic boundary conditions
  PBC(Pos)
  allvar[i]=DiffVar(Pos,InitPos) #calculate diffusion at current step
# recalculate forces
  fattdt=Forces(Pos)
# get velocities
  for j in xrange(Npart):
   VelXYZ[j][0]=VelXYZ[j][0]+0.5*(fattdt[j][0]+fatt[j][0])*dTiSt # x-component
   VelXYZ[j][1]=VelXYZ[j][1]+0.5*(fattdt[j][1]+fatt[j][1])*dTiSt # y-component
   VelXYZ[j][2]=VelXYZ[j][2]+0.5*(fattdt[j][2]+fatt[j][2])*dTiSt # z-component
  U,sfxr,rijs=Epot(Pos)
  allFixRi[i]=sfxr
  T=Ekin(VelXYZ)
  KinEgs[i]=T
  temp=T/1.5/Npart
  print 'time step %5d temperature= %5.2f %s Etot= %12.7f %s' %(i+1,temp*eps, 'K', U+T,'in MD units')
  cormom(VelXYZ,Npart)
  scf= mt.sqrt(T/Ekin(VelXYZ)) #scaling factor
  VelXYZ=scf*VelXYZ 
  Save(Pos,TiStPos,i+1,rijs,allvar) #save positions & generate the radial distribution function 
 return KinEgs, allFixRi

# Nose-Hoover thermostat
def ThermoNH(Pos,VelXYZ,dTiSt,tempD):
 temp=Ekin(VelXYZ)/1.5/Npart
 beta=mt.sqrt(tempD/temp)
 zeta=0. # initial guess for zeta (friction parameter)
 thresNH=0.00001 # threshold in kelvins
 MNH=0.00001 # coupling constant for the Nose-Hoover thermostat 
 jCycle=0
 print 'adjusting temperature to %6.3f K, expected deviation is %6.3f K' %(eps*tempD, mt.sqrt(2./3./Npart)*eps*tempD)
 while eps*abs(temp-tempD)>thresNH:
  jCycle=jCycle+1 
  fatt=Forces(Pos)
  for j in xrange(Npart):
   pref0=1-0.5*dTiSt*zeta 
# update positions   
   Pos[j][0]=Pos[j][0]+pref0*VelXYZ[j][0]*dTiSt+0.5*fatt[j][0]*dTiSt**2 # x-component
   Pos[j][1]=Pos[j][1]+pref0*VelXYZ[j][1]*dTiSt+0.5*fatt[j][1]*dTiSt**2 # y-component 
   Pos[j][2]=Pos[j][2]+pref0*VelXYZ[j][2]*dTiSt+0.5*fatt[j][2]*dTiSt**2 # z-component
# periodic boundary conditions
  PBC(Pos)
# update zeta
  oldzeta=zeta
  zeta=zeta+dTiSt*(2*Ekin(VelXYZ)-3*Npart*tempD)/MNH
  pref1=(1-0.5*dTiSt*oldzeta)/(1+0.5*dTiSt*zeta) # prefactor 1 for updating velocities
  pref2=1./(1+0.5*dTiSt*zeta) # prefactor 2 for updating velocities
# recalculate forces
  fattdt=Forces(Pos)
# get updated velocities
  for j in xrange(Npart):
   VelXYZ[j][0]=pref1*VelXYZ[j][0]+0.5*pref2*(fattdt[j][0]+fatt[j][0])*dTiSt # x-component
   VelXYZ[j][1]=pref1*VelXYZ[j][1]+0.5*pref2*(fattdt[j][1]+fatt[j][1])*dTiSt # y-component
   VelXYZ[j][2]=pref1*VelXYZ[j][2]+0.5*pref2*(fattdt[j][2]+fatt[j][2])*dTiSt # z-component
  cormom(VelXYZ,Npart)
  KE=Ekin(VelXYZ)
  scf= mt.sqrt(KE/Ekin(VelXYZ)) #scaling factor
  VelXYZ=scf*VelXYZ 
  temp=KE/1.5/Npart
  print 'step %5d current temperature = %6.3f K dtemp = %6.3f K zeta = %6.3f oldzeta = %6.3f' %(jCycle,temp*eps,eps*(tempD-temp),zeta,oldzeta)
 return

# Nose-Hoover thermostat 2
def ThermoNH2(Pos,VelXYZ,dTiSt,tempD):
 temp=Ekin(VelXYZ)/1.5/Npart
 beta=mt.sqrt(tempD/temp)
 zeta=0. # initial guess for zeta (friction parameter)
 thresNH=0.1 # threshold in kelvins
 MNH=1 # coupling constant for the Nose-Hoover thermostat 
 jCycle=0
 VelXYZ12=np.empty(3*Npart).reshape(Npart,3)
 print 'adjusting temperature to %6.3f K, expected deviation is %6.3f K' %(eps*tempD, mt.sqrt(2./3./Npart)*eps*tempD)
 while eps*abs(temp-tempD)>thresNH:
  jCycle=jCycle+1
  fatt=Forces(Pos) 
  for j in xrange(Npart):
   pref0=1-0.5*dTiSt*zeta 
# update positions   
   Pos[j][0]=Pos[j][0]+pref0*VelXYZ[j][0]*dTiSt+0.5*fatt[j][0]*dTiSt**2 # x-component
   Pos[j][1]=Pos[j][1]+pref0*VelXYZ[j][1]*dTiSt+0.5*fatt[j][1]*dTiSt**2 # y-component 
   Pos[j][2]=Pos[j][2]+pref0*VelXYZ[j][2]*dTiSt+0.5*fatt[j][2]*dTiSt**2 # z-component
# periodic boundary conditions
  PBC(Pos)
  oldzeta=zeta
  for j in xrange(Npart):
   VelXYZ12[j][0]=VelXYZ[j][0]+0.5*dTiSt*(fatt[j][0]-oldzeta*VelXYZ[j][0]) # x-component
   VelXYZ12[j][1]=VelXYZ[j][1]+0.5*dTiSt*(fatt[j][1]-oldzeta*VelXYZ[j][1]) # x-component
   VelXYZ12[j][2]=VelXYZ[j][2]+0.5*dTiSt*(fatt[j][2]-oldzeta*VelXYZ[j][2]) # x-component
# update zeta
  zeta=zeta+dTiSt*(2*Ekin(VelXYZ12)-3*Npart*tempD)/MNH
# recalculate forces
  fattdt=Forces(Pos)
# get updated velocities
  pref=1./(1+0.5*dTiSt*zeta)
  for j in xrange(Npart):
   VelXYZ[j][0]=pref*(VelXYZ12[j][0]+0.5*dTiSt*fattdt[j][0])# x-component
   VelXYZ[j][1]=pref*(VelXYZ12[j][1]+0.5*dTiSt*fattdt[j][1])# y-component
   VelXYZ[j][2]=pref*(VelXYZ12[j][2]+0.5*dTiSt*fattdt[j][2])# y-component
#  cormom(VelXYZ,Npart)
  KE=Ekin(VelXYZ)
#  scf= mt.sqrt(KE/Ekin(VelXYZ)) #scaling factor
#  VelXYZ=scf*VelXYZ 
  temp=KE/1.5/Npart
  print 'step %5d current temperature = %6.3f K dtemp = %6.3f K zeta = %6.3f oldzeta = %6.3f' %(jCycle,temp*eps,eps*(tempD-temp),zeta,oldzeta)
 return

# print summary of the calculation
def summary(KEs,sumFR,Npart,dTiSt,NTiSt,Lbox):
 KE=np.average(KEs) # average kinetic energy
 sumFR=np.average(sumFR) 
 T=KE/1.5/Npart #average temperature
 P=1/(3.*Lbox**3)*(2*KE+1./3.*sumFR)
 gbyemsg="""
 --------------------------------------------------------------------------------
 SUMMARY\n
 total simulation time = %5.2f fs
 average temperature = %5.2f K
 average kinetic energy = %10.5f in MD units
 pressure = %5.5f atm
 --------------------------------------------------------------------------------
 """ %(dTiSt*NTiSt*2.068*10**3,T*eps,1.5*Npart*T,P*eps*sp.constants.k/((sigma*10**(-10))**3)/101325.)
 print gbyemsg
 return

# calculate diffusion (mean standard deviation of the particle positions)
def DiffVar(Pos,InitPos):
 sumvars=0.
 for i in xrange(Npart):
  vari=(sp.spatial.distance.euclidean(Pos[i],InitPos[i]))**2
  sumvars=sumvars+vari
 return sumvars

# generate radial distibution functions
def RadDistr(rijs):
 data=np.array([],dtype=float).reshape(0,2)
 dr=0.001 # step for histogram
 nsteps=mt.ceil(np.max(rijs)/dr)
 rad=np.linspace(0,max(rijs),num=nsteps) #intervals
 nint=len(rad)-1
 for i in xrange(nint):
  dv=4./3.*mt.pi*(rad[i+1]**3-rad[i]**3) 
  data=np.vstack([data, [rad[i+1], len(rijs[(rijs>=rad[i]) & (rijs<=rad[i+1])])/dv]])
 return nint,data

# dump coordinates of particles at certain time steps
def Save(Pos,TiStPos,TiSt,rijs,allvar):
# if len(TiStPos) == 0:
#  return
 for j in TiStPos:
  if TiSt == j:
   j=str(j)
# positions
   filename=j+'.pos' # file to store positions of particles at jth step
   outfile=open(filename,'w')
   header="""#INPUT DATA
#
#Parameters of the system: 
#A number of particles = %5d
#Dimension of the box = %4.2f (m*10**10) = %4.2f (m.d.u)
#Mass of a particle = %6.3f (a.u.)
#Temperature = %5.2f K
#
#Integration (velocity-Verlet):
#Time step = %6.3f fs
#A number of time steps = %5d
#
#Lennard-Jones (LJ) potential:
#Epsilon = %4.1f K (in units of k_B)
#Sigma = %6.3f (m*10**10)
#Cutoff radius for interparticle interactions = %4.2f (in units of the LJ sigma)
#Switching parameter for the potential = %4.2f (in units of the LJ sigma)
#
#Coordinates of particles at time step %4s
#X Y Z\n""" % (Npart,Lbox,Lbox/sigma,mass,T,dTiSt,NTiSt,eps,sigma,Rcut,Rl,j)
   outfile.write(header)  
   for k in xrange(Npart):
    outfile.write( '%12.6f %12.6f %12.6f\n' %(Pos[k][0],Pos[k][1],Pos[k][2])) 
   outfile.close()
# pair-correlation function
   nint,gr=RadDistr(rijs)
   filename=j+'.rd' # file to store pair-correlation function at jth step
   outfile=open(filename,'w')
   header="""#INPUT DATA
#
#Parameters of the system: 
#A number of particles = %5d
#Dimension of the box = %4.2f (m*10**10) = %4.2f (m.d.u)
#Mass of a particle = %6.3f (a.u.)
#Temperature = %5.2f K
#
#Integration (velocity-Verlet):
#Time step = %6.3f fs
#A number of time steps = %5d
#
#Lennard-Jones (LJ) potential:
#Epsilon = %4.1f K (in units of k_B)
#Sigma = %6.3f (m*10**10)
#Cutoff radius for interparticle interactions = %4.2f (in units of the LJ sigma)
#Switching parameter for the potential = %4.2f (in units of the LJ sigma)
#
#Pair-correlation function at time step %4s
#R G(R)\n""" % (Npart,Lbox,Lbox/sigma,mass,T,dTiSt,NTiSt,eps,sigma,Rcut,Rl,j)
   outfile.write(header)  
   for k in xrange(nint):
    outfile.write( '%6.3f %12.6f\n' %(gr[k][0],gr[k][1])) 
   outfile.close()
# diffusion
 if TiSt == NTiSt:
  filename='diff' # file to store pair-correlation function at jth step
  outfile=open(filename,'w')
  header="""#INPUT DATA
#
#Parameters of the system: 
#A number of particles = %5d
#Dimension of the box = %4.2f (m*10**10) = %4.2f (m.d.u)
#Mass of a particle = %6.3f (a.u.)
#Temperature = %5.2f K
#
#Integration (velocity-Verlet):
#Time step = %6.3f fs
#A number of time steps = %5d
#
#Lennard-Jones (LJ) potential:
#Epsilon = %4.1f K (in units of k_B)
#Sigma = %6.3f (m*10**10)
#Cutoff radius for interparticle interactions = %4.2f (in units of the LJ sigma)
#Switching parameter for the potential = %4.2f (in units of the LJ sigma)
#
#Diffusion (mean standard deviation of the particle positions)
#TiSt Var(TiSt)\n""" % (Npart,Lbox,Lbox/sigma,mass,T,dTiSt,NTiSt,eps,sigma,Rcut,Rl)
  outfile.write(header)  
  for k in xrange(NTiSt):
   outfile.write( '%4i %20.10f\n' %(k+1,allvar[k])) 
  outfile.close()
 return

#----------------------------------------------------------------------------
# where the code starts...
# initial parameters
Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt,TiStPos,Rcut,Rl=input() 
iout(Npart,Lbox,mass,eps,sigma,T,dTiSt,NTiSt)

# conversion
Lbox,T,dTiSt=convr(Lbox,T,dTiSt,mass,eps)

# generate the grid
Pos=GenGrd(Npart,Lbox)
InitPos=np.empty(Npart*3).reshape(Npart,3) #needed for calculation of diffusion

Ek=1.5*Npart*T

print '%s Etot=%12.7f %s' %('initial energy:',Epot(Pos)[0]+Ek,'in MD units')

# generate initial velocities
VelXYZ=GenVel(T)

#ThermoNH(Pos,VelXYZ,dTiSt,4*T)
ThermoNH2(Pos,VelXYZ,dTiSt,100.)

# integrator 
#KEs,sumFR=VelVer(Pos,InitPos,VelXYZ,dTiSt,NTiSt)
#summary(KEs,sumFR,Npart,dTiSt,NTiSt,Lbox)
