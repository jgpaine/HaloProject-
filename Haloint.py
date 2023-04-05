import numpy as np
from scipy import integrate

def rho(r):
    ps=1.0
    rs=4.0
    gamma=3 
    #return ps*(rs/r)*(1+r/rs)**-2 #NFW 
    return 1.453*np.exp(-r**0.45)
    #return 2**(3-3*gamma)*ps * 1/ ((rs/r)^gamma) * 1/((1+r/rs)**-2)


def integrand( x, b):
    return rho(np.sqrt(x**2 + b**2) ) * b 

def compute_integral(rho_func, b, x, bil):
    integrand_func = np.vectorize(lambda x, b: integrand(x,b))
    integral, _ = integrate.dblquad(integrand_func, 0, bil, -np.inf, np.inf)
    return 2 * np.pi * integral

def deflection(Mass, Dl, Di, bil):
    G= 000.43
 #  4.3009172706(3)×10−3	pc⋅M⊙−1⋅(km/s)2
    deflection = -(1- (Dl/Di))*G*Mass/bil
    return deflection


# Example usage, testing 
b = np.array( [ 3.0])
x = np.array([ 3.0])
bil = 5.0


x_array = np.array([ 1,4,5,6,7,8])  
b_array = np.array([ 1,4,5,6,7,8])  
EnclosedMass = np.empty((len(x_array), len(b_array)))
DeflectionAngles= np.empty((len(x_array), len(b_array)))

for i in range(1,len(x_array)):
    for j in range(1,len(b_array)):
        integral = compute_integral(rho, b_array[j], x_array[i], bil)
        EnclosedMass[i, j] = integral
        DeflectionAngles[i,j] = deflection( EnclosedMass[i, j], 10, 1, bil)



print("x", "\t", "b", "\t", "Enclosed Mass", " \t", "Angle of deflection")
for i, xi in enumerate(x_array):
    for j, bj in enumerate(b_array):
        print(xi, "\t", bj, "\t", EnclosedMass[i, j],DeflectionAngles[i,j] ) 



