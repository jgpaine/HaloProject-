import numpy as np
from scipy import integrate
from scipy.integrate import tplquad

## This works !! 

def rho(r): # This funcation is the mass density profile of the dark matter subhalo
    ps=1.0
    rs=4.0
    gamma=3 
    #return ps*(rs/r)*(1+r/rs)**-2 #NFW 
    return 1.453*np.exp(-r**0.45)
    #return 2**(3-3*gamma)*ps * 1/ ((rs/r)^gamma) * 1/((1+r/rs)**-2)

def integrand( x, b): #This is the integrand to get the enclosed mass from Halometry Paper Tilburg
    return rho(np.sqrt(x**2 + b**2) ) * b 

def compute_integral(rho_func, b, x): # this integral get the enclosed mass at point (x,b)
    bil = np.sqrt(x**2 + b**2) 
    integrand_func = np.vectorize(lambda x, b: integrand(x,b))
    integral, _ = integrate.dblquad(integrand_func, 0, bil, -np.inf, np.inf)
    return 2 * np.pi * integral

def deflection(Mass, Dl, Di, bil):
 # calulates the deflection angle of a photon due the mass enclosedm,the distance to lens, disntance to images/source, and the impact paremter of the photon bil
    G= 000.43
    # galaxy 20 kpc
 #  4.3009172706(3)×10−3	pc⋅M⊙−1⋅(km/s)2
    deflection = -(1- (Dl/Di))*G*Mass/bil # From  equation 2.1 in Halometery
    #want in mico arc seconds
    #check with 2.70 and 2.39 in dodelson
    return deflection
    
## add graviational potential intergral
def gravitational_potential(rho_func, r_vals):
    """
    Calculates the gravitational potential from a mass density function rho(r).
    
    Parameters:
        rho_func (function): A function that takes in a radial coordinate r and returns the mass density rho(r).
        r_vals (array): An array of radial coordinates at which to evaluate the potential.
        
    Returns:
        phi_vals (array): An array of the gravitational potentials evaluated at the radial coordinates r_vals.
    """
    G = 000.43 # Gravitational constant
    phi_vals = np.zeros_like(r_vals)  # Initialize the potential array
    for i, r in enumerate(r_vals):
        integrand = lambda x: G * rho_func(x) * x**2  / np.sqrt(x**2 - r**2)
        integral, _ = integrate.quad(integrand, r, np.inf)
        phi_vals[i] =  integral
    return phi_vals
    
#print(gravitational_potential(rho, [4,5,6,7,7]))


def Magofpoint (B,theta,x,b):
# 4.16 for a point mass, B is eq.5.2 need to inetrgrate over line of cite to get
     bil = np.sqrt(x**2 + b**2) 
     # check what beta is in text 
     u = B**2 + 2 * theta**2 / (B* np.sqrt(B**2 + 4*theta**2))
     # book says angle for enstien ring, but I think will change with the change in angle
     # not sure if this should be a constant or not
     return u

def Boverlineofsight(Di,z,Phi): #Eq. 2.52
#Phi to be an array of graviatonal pontentails
# equation used deriavtive of graviational pontential, think I need this for magnification generally
# yay don' t need this!!
  return B

def Mag(k,gamma1,gamma2)
# K is convergernce and shear is gamma1 and gamma1
    u= [(1-k**2)-(gamma1**2+gamma2**2))]**-1
    return u

Def Kappa(phi,thetax,thetay):
    
 return k

## need general equation for magificantion
# Example usage, testing 
#b = np.array( [ 3.0])
#x = np.array([ 3.0])
#bil = 5.0


"""
This now runs the funcations we made. In putting an array of x and b locations of a photon
and using the density profile rho, gets the mass enclosed and angle of deflection at many points
"""
x_array = np.array([ 1,5,10,30,50,100])
b_array = np.array([ 1,4,5,6,7,8])  
EnclosedMass = np.empty((len(x_array), len(b_array)))
DeflectionAngles= np.empty((len(x_array), len(b_array)))

for i in range(1,len(x_array)):
    for j in range(1,len(b_array)):
        integral = compute_integral(rho, b_array[j], x_array[i])
        EnclosedMass[i, j] = integral
        bil = np.sqrt(x_array[i]**2 + b_array[j]**2)
        # in put distnace dis to len, dis to image 
        DeflectionAngles[i,j] = deflection( EnclosedMass[i, j], 15, 5, bil)



print("x", "\t", "b", "\t", "Enclosed Mass", " \t", "Angle of deflection")
for i, xi in enumerate(x_array):
    for j, bj in enumerate(b_array):
        print(xi, "\t", bj, "\t", EnclosedMass[i, j],DeflectionAngles[i,j] ) 



