#!/usr/bin/env python3


import numpy as np
from   scipy.optimize import leastsq
import scipy.constants as sc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import math


Bohr2Ang = sc.physical_constants['Bohr radius'][0]*1e10
Ry2eV = sc.physical_constants['Rydberg constant times hc in eV'][0] 
qe = sc.physical_constants['elementary charge'][0]


#=====================================


def emto_read_post_data():

    filename='emto_post_data'
    

    jobn=[]
    Etot=[]
    SWS=[]

    f = open(filename, 'r')
    next(f)
    for line in f:
        if np.abs( float(line.split( )[1]) ) > 1e-6 : 
            jobn.append( float(line.split( )[0]) )
            Etot.append( float(line.split( )[1]) )
            SWS.append(  float(line.split( )[2]) )
    
    print('\nCHECK jobn Etot(Ry) SWS(Bohr):')
    print(jobn)
    print(Etot)
    print(SWS)
    
    f.close()

    return jobn, Etot, SWS


jobn, Etot, SWS=emto_read_post_data()



if len(SWS) < 3.9 :
    import sys
    sys.exit("\nABORT: less than 4 data points!\n")



#=====================================


vols = []
energies = []

for i in range(0, len(SWS), 1):
    energies.append( Etot[i]*Ry2eV )
    vols.append( 4/3*math.pi*(SWS[i]*Bohr2Ang)**3 )


jobn=np.array(jobn)
vols=np.array(vols)
energies=np.array(energies)


print('\nCHECK jobn energies(eV) vols(Ang^3):')
print(jobn)
print(energies)
print(vols)


def myeos(parameters, vol):
    E0, B0, BP, V0 = parameters

# Birch-Murnaghan equation of state
    E= E0 + 9*V0*B0/16 *( \
    (  (V0/vol)**(2/3)-1 )**3 *BP \
    +( (V0/vol)**(2/3)-1 )**2 *( 6- 4*(V0/vol)**(2/3) ) \
    )

# From Phys. Rev. B 28, 5480 (1983)
#    E = E0 \
#    + B0 * vol / BP * (((V0 / vol)**BP) / (BP - 1) + 1) \
#    - V0 * B0 / (BP - 1.0)

    return E



# we will minimize this function
def myerrfunc(pars, y, x):
    err =  y - myeos(pars, x)
    return err


# initial guess of parameters
x0 = [ min(energies), 1.0, 5.0, np.mean(vols)]

plsq, cov, infodict, mesg, ier = \
leastsq(myerrfunc, x0, args=(energies, vols), full_output=1 )

print( '\nFitted parameters = {0}\n'.format(plsq) )



ssErr = (infodict['fvec']**2).sum()
ssTot = ((energies-energies.mean())**2).sum()
R2 = 1-(ssErr/ssTot )



E0=plsq[0]
B0=plsq[1]*qe*1e21
BP=plsq[2]
V0=plsq[3]

SWS0=(V0/(4/3*math.pi))**(1/3) /Bohr2Ang 

a0=(V0*4)**(1/3) 

V1V0= np.mean( np.divide( vols/V0, jobn ))

#=====================================


f = open("emto_Birch_Murnaghan_EOS.txt","w+")

f.write("# EOS fitting: \n" )

f.write("\n%16s %16s %16s %16s \n" \
%('E0 (eV)', 'V0 (Ang^3)', 'B0 (GPa)', 'B1') )

f.write("%16.8f %16.8f %16.8f %16.8f \n" \
%(E0, V0, B0, BP) )


f.write("\n%16s %16s %16s %16s \n" \
%('SWS0 (Bohr)', 'a0_fcc (Ang)', 'R2', 'V1/V0-1') )

f.write("%16.6f %16.8f %16.8f %16.8f \n" \
%( SWS0, a0, R2, (V1V0-1) ) )


f.write("\n%16s %16s \n" %('E (eV)', 'V (Ang^3)') )
for i in range(0, len(vols), 1):
     f.write("%16.8f %16.8f \n" %( energies[i], vols[i]) )

f.write("\n%16s %16s \n" %('E (Ry)', 'SWS (Bohr)') )
for i in range(0, len(vols), 1):
     f.write("%16.8f %16.8f \n" %( Etot[i], SWS[i]) )

f.close() 


#=====================================


def cm2inch(value):
    return value/2.54


plt.rcParams['font.size']=10
#plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth']=0.5
plt.rcParams['axes.grid']=True
plt.rcParams['grid.linestyle']='--'
plt.rcParams['grid.linewidth']=0.5
plt.rcParams["savefig.transparent"]='True'

plt.rcParams['lines.linewidth']=0.8
plt.rcParams['lines.markersize'] = 5


fig = plt.figure(figsize=( cm2inch(8), cm2inch(7) ))
ax = fig.add_axes([0.18, 0.18, 0.78, 0.78 ])


plt.plot(vols/V0, (energies-E0)*1e3, '+', \
color=np.array([ 223,   0,  36])/255 )


#plot the fitted curve on top
x = np.arange(min(vols), max(vols), 0.0001*V0)
y = myeos(plsq, x)
plt.plot(x/V0, (y-E0)*1e3, '-', color=np.array([ 0, 159,  61])/255  ) 


plt.xlabel('$V/V_0$')
plt.ylabel('$E-E_0$ (meV)')
plt.xticks(np.arange(0.8, 1.2, 0.05))
plt.xlim([0.92, 1.08])
plt.savefig('emto_Birch_Murnaghan_EOS.pdf')


