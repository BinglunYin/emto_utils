#!/usr/local/bin/python3

from sympy import *
import numpy as np 

def main():
    e0, d= symbols('e0 d')

    # pre-strain
    D0 = Matrix([[1+e0, 0, 0], [0, 1+e0, 0], [0, 0, 1+e0]])
   
    # orth, mono
    Do = Matrix([[1+d, 0, 0], [0, 1-d, 0], [0, 0, 1/(1-d**2)]])
    Dm = Matrix([[1, d, 0], [d, 1, 0], [0, 0, 1/(1-d**2)]])
    
    # Do = Matrix([[1+d, 0, 0], [0, 1, 0], [0, 0, 1]])
    print('\n==> det:')
    print(simplify( D0.det() ))
    print(simplify( Do.det() ))
    print(simplify( Dm.det() ))

    # orth
    eo = simplify( D0 * Do - diag(1, 1, 1) )
    Eo = calc_E_from_strain(eo)

    print('\n==> with e0=0, Eo=:')
    print(Eo.subs([(e0, 0)]))

    
    # mono
    em = simplify( D0 * Dm - diag(1, 1, 1) )
    Em = calc_E_from_strain(em)

    print('\n==> with e0=0, Em=:')
    print(Em.subs([(e0, 0)]))





def calc_E_from_strain(e):
    # for cubic 
    C11, C12, C44, d, e0 = symbols('C11 C12 C44 d e0')

    E = 1/2*C11 *( e[0,0]**2 +e[1,1]**2 +e[2,2]**2 ) \
           +C12 *( e[0,0]*e[1,1] +e[0,0]*e[2,2] +e[1,1]*e[2,2] )  \
       +  2*C44 *( e[0,1]**2 +e[0,2]**2 +e[1,2]**2 )  
    
    E = expand(E)
    E = collect(E, d)
      
    print('\n==> energy:'  )
    print(E)
    for i in np.arange(5):
        print('\n==> coefficient of d** %1d' %(i)  )
        print(collect( E.coeff(d**i), e0) )

    d0 = simplify( -E.coeff(d) /E.coeff(d**2) /2 ) 
    print('\n==> -b/2a='  )
    print(collect(d0, e0))
    
    return E


main()



