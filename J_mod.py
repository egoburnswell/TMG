"""
GF_mod contains a range of functions related to the Green's functions of graphene and nanotubes and the coupling calculations between impurities.

GLOBAL VARIABLES:
t0 is the unstrained Graphene or Nanotube hopping integral.
tau is the hopping integral between the impurity and the substrate.
t1, t2 are the strained hopping integrals for Graphene and Nanotubes.
Ef is the Fermi Energy.
eta is the small imaginary part added to calculate the Green's Function.
i is used in place of the standard python imaginary number 1j.

OTHER VARIABLES:
m is a multiplier to the vector a1.
n is a multiplier to the vector a2.
S1 takes an integer value 1 or 2 (or 3 or 4 in the 4-atom unit cell), and corresponds to the sublattice of site A.
S2 takes an integer value 1 or 2 (or 3 or 4 in the 4-atom unit cell), and corresponds to the sublattice of site B.
En is an abbreviation for Energy.
Nc is the width of a nanotube.
"""

# IMPORTS
import GF_MOD
from scipy.integrate import quad
from numpy           import matrix, zeros, transpose, linspace
from numpy.linalg    import inv, det
from numpy.matlib    import identity
from cmath           import sin, cos, log, acos, sqrt, exp, pi
from time            import time
from numpy		import random
from functools import partial
import multiprocessing

# VARIABLES
global  t0, t1, t2, tau, Ef, i, eta
i=1j
t0 = -1.0
t1 = -1.0
t2 = -1.0
eta = 0.0001

# FORTRAN RELATED FUNCTIONS 
def C_int(f,lim1,lim2, eps, lmt):		# This is nicer than your current python one. Keep it.  
  int_re = quad(lambda x: f(x).real, lim1, lim2, epsabs=0.0, epsrel=eps, limit=lmt )
  int_im = quad(lambda x: f(x).imag, lim1, lim2, epsabs=0.0, epsrel=eps, limit=lmt )
  return int_re[0] + 1j*int_im[0]


def Green(Gf_eps, Gf_lmt,E, m,n,s_lat):
  """The Graphene Green's function
    The kZ integration is performed last"""
  GF = partial(GF_MOD.gf_bulk_kz,m,n,s_lat,E)
  return C_int(GF,-pi/2,0, Gf_eps, Gf_lmt) + C_int(GF,0,pi/2, Gf_eps, Gf_lmt)
  #return C_int(GF,-pi/2, pi/2, Gf_eps, Gf_lmt) 

def Green_strain(Gf_eps, Gf_lmt,E, m, n, s ):
  """The Graphene Green's function
    The kZ integration is performed last"""
  GF = partial( GF_MOD.gf_bulk_kz_strain, t1, t2, m, n, s, E )
  return C_int(GF,-pi/2,0, Gf_eps, Gf_lmt) + C_int(GF,0,pi/2, Gf_eps, Gf_lmt)
# FORTRAN RELATED FUNCTIONS 


def set_strain( strain = 0, angle = 0 ):
    """  
    Sets Global values to be used in all calculations.
    Strain = 0.0 - 0.2 for realistic strain.
    Angle: 0 = Zigzag direction, pi/2 = Armchair direction, only works for these directions!!
    """  
    global t1, t2
    # Strain matrix
    sig = 0.165
    E11 = strain*(cos(angle)**2 - sig*sin(angle)**2)
    E12 = strain*((1+sig)*cos(angle)*sin(angle))
    E22 = strain*(sin(angle)**2 - sig*cos(angle)**2)
    # New distances
    l1 = 1 + E22
    l2 = 1 + 0.75*E11 - (sqrt(3.0)/2)*E12 + 0.25*E22
    # New hoppings
    t1 = -exp(-3.37*(l1-1)) 
    t2 = -exp(-3.37*(l2-1))
    
    




"""For use with testing code. """
if __name__ == "__main__":
  from numpy import arange
  
## Green's function test  
  set_strain( strain = 0.1 )
  for E in arange(-3,3, 0.01):
     En = E + i*eta
     G2 = Green_strain( 1.0e-4, 100,  En, 0, 0, 0)
     print En.real, G2.real, G2.imag
      



