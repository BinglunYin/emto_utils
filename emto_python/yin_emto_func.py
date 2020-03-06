#!/usr/local/bin/python3


import numpy as np



def phy_const(sym):
    pc={}
    pc['Bohr2Ang'] = 5.29177210903e-1      
    pc['Ry2eV'] = 13.605693112994   
    pc['qe'] = 1.602176634e-19 

    y = pc[sym]
    return y



def emto_read_post_data(filename='emto_post_data'):
    jobn=np.array([])
    Etot=np.array([])
    SWS =np.array([])

    f = open(filename, 'r')
    next(f)
    for line in f:
        # to skip unfinished job
        if np.abs( float(line.split( )[1]) ) > 1e-6 :  
            jobn = np.append( jobn, float(line.split()[0]) )
            Etot = np.append( Etot, float(line.split()[1]) )
            SWS  = np.append( SWS,  float(line.split()[2]) )
    f.close()

    print('==> CHECK jobn Etot(Ry) SWS(Bohr):')
    print(jobn, Etot, SWS)
    return jobn, Etot, SWS






