#!/home/yin/opt/bin/python3


import numpy as np
import yin_emto_func as ef


def main():
    jobn1, Etot1, SWS1 = ef.emto_read_post_data()
    x1, y1 = post_data( jobn1, Etot1, SWS1 )
    fitres1 = myfitting(x1, y1)

    jobn2, Etot2, SWS2 = ef.emto_read_post_data('../emto_Cij_mono/emto_post_data')
    x2, y2 = post_data( jobn2, Etot2, SWS2 )
    fitres2 = myfitting(x2, y2)
    
    B = read_B()

    Cij = np.array([
        (2*fitres1[0] +3*B)/3,
        ( -fitres1[0] +3*B)/3,
           fitres2[0]/2
    ])
    print(Cij)

    write_output(Cij, fitres1, fitres2)
    plot_Cij(x1, y1, x2, y2, fitres1, fitres2)




#=====================================


def post_data(jobn, Etot, SWS):
    n=2

    if len(jobn) < n-0.1:
        import sys
        sys.exit("\n==> ABORT: less than %d data points!\n" %(n))
 
    Bohr2Ang = ef.phy_const('Bohr2Ang')  
    Ry2J = ef.phy_const('Ry2eV') * ef.phy_const('qe') 

    V = 4/3 *np.pi *(SWS*Bohr2Ang)**3 
    V0 = V.mean()
  
    x = jobn.copy()            #[strain]
    y = Etot * Ry2J /V0 *1e21  #[GPa]
    print(x, y)
    return x, y



#==========================

# fitting eqn
def myeqn(param, x):
    a, c  = param
    y= a* x**2 +c 
    return y 
    


# minimize this function
def myerrfunc(param, y, x):
    err =  y - myeqn(param, x)
    return err

   

def myfitting(x, y):
    param0 = [ 100, -6e5]
   
    from scipy.optimize import leastsq

    plsq, cov, infodict, mesg, ier = \
    leastsq(myerrfunc, param0, args=(y, x), full_output=1 )
   
    print( '\n==> fitted parameters = {0}\n'.format(plsq) )

    ssErr = (infodict['fvec']**2).sum()
    ssTot = ((y-y.mean())**2).sum()
    R2 = 1-(ssErr/ssTot )

    fitres=np.array([plsq[0], plsq[1], R2])

    return fitres



#=====================================

def read_B():
    filename = '../emto_E_V/emto_post_eos.txt'
    f = open(filename, 'r')
    next(f);  next(f)
    line1 = f.readline()
    B = float(line1.split()[2]) 
    f.close()
    return B



#=====================================

def write_output(Cij, fitres1, fitres2):

    f = open("emto_post_Cij_cubic.txt","w+")

    f.write("# emto Cij cubic: \n" )

    f.write("%16s %16s %16s \n" \
    %('C11 (GPa)', 'C12 (GPa)', 'C44 (GPa)') )

    f.write("%16.8f %16.8f %16.8f \n" \
    %(Cij[0], Cij[1], Cij[2]) )

    f.write("\n%16s %16s \n" \
    %('(C11+2C12)/3', '(C11-C12)/2' ) )

    f.write("%16.8f %16.8f \n" \
    %( (Cij[0] + 2*Cij[1])/3 , (Cij[0]-Cij[1])/2      ) )


    f.write("\n# fitting results of y = a*x**2 +c [GPa]: \n" )

    f.write("\n%16s %16s %16s \n" \
    %('a1 (orth)', 'c1', 'R21') )

    f.write("%16.8f %16.8e %16.8f \n" \
    %( fitres1[0], fitres1[1], fitres1[2] ) )   

    f.write("\n%16s %16s %16s \n" \
    %('a2 (mono)', 'c2', 'R22') )

    f.write("%16.8f %16.8e %16.8f \n" \
    %( fitres2[0], fitres2[1], fitres2[2] ) )   

    f.close() 



#=====================================

def plot_Cij(x1, y1, x2, y2, fitres1, fitres2):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams['font.size']=8
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.linewidth']=0.5
    plt.rcParams['axes.grid']=True
    plt.rcParams['grid.linestyle']='--'
    plt.rcParams['grid.linewidth']=0.2
    plt.rcParams["savefig.transparent"]='True'
    plt.rcParams['lines.linewidth']=0.8
    plt.rcParams['lines.markersize'] = 4
 

    fig_w = 3.15
    fig_h = 5

    xi=np.arange(0, 0.0501, 1e-4)

    fig1, ax1 = plt.subplots(nrows=2, ncols=1, \
    sharex=True, figsize=(fig_w, fig_h) )

    ax1[0].plot(x1, (y1-fitres1[1])*1e3, 'o')
    yi1 = fitres1[0] *xi**2 *1e3
    ax1[0].plot(xi, yi1, '-')

    ax1[1].plot(x2, (y2-fitres2[1])*1e3, 'o')
    yi2 = fitres2[0] *xi**2 *1e3
    ax1[1].plot(xi, yi2, '-')

    pos=np.array([0.18, 0.55, 0.75, 0.4])
    ax1[0].set_position(pos)
    ax1[1].set_position(pos +np.array([0, -0.45, 0, 0]) )

    plt.setp(ax1[-1], xlabel='$\delta$')
    plt.setp(ax1[:],  ylabel='energy density (MPa)')

    ax1[0].text(0.0005, round(yi1.max()*0.8), \
    'orth: a= %.1f GPa'%(fitres1[0])  )

    ax1[1].text(0.0005, round(yi2.max()*0.8), \
    'mono: a= %.1f GPa'%(fitres2[0])  )

    plt.savefig('emto_post_Cij_cubic.pdf')




main()





