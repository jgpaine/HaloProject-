import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# june 6 2023, Make the plot in the book


mass = 10*10**25  # solar mass
# Distances in kpc around 10 kps 1 kiloparsec = 3.086 * 10**16
Dl = 1 * 3.086 * 10**16
Dsl = 1 * 3.086 * 10**16
Ds = Dl + Dsl
v = 100 # km/s
c2= 9*10**10 #km/s^2 
# Gravitational constant in km^2 Mpc M_sun^-1 s^-2
G = 4.301 * 10**-9
# add converstion

# B values in km
bx = [-9,-8,-7,-6,-5,-4,-3,-2,-1,-0.5,-0.25,-0.05,0.000001,0.05,0.25,0.5,1, 2, 3, 4, 5, 6, 7, 8, 9]
# define this in mass loop
def thetaE(Mass, Grav, Dstol, Dtos, Dtol): # eq.2.70
 
    thetaE = math.sqrt((4 * Grav * Mass * Dstol) / (c2 * Dtos * Dtol * Dstol)) # add c^2
    return thetaE

def total_magnification(Beta, thetaE): # eq 4.16
    u = Beta**2 + 2 * thetaE**2 / (Beta * np.sqrt(Beta**2 + 4 * thetaE))
    return u

def time(D,th,velocity): # 5.4  use 5.5
    t =( (D*th) / velocity ) 
    return t
def duration(Mass,DL,DSL,DS,v):
    td= 78 * (200/v)* Mass**0.5 *(DL/10)**0.5 * (DSL/DS)**0.5
    return td

def compute_mu(bx,by, theta):
    u = (math.sqrt(bx**2 + by**2)) / theta
    numerator = u**2 + 2
    denominator = u * math.sqrt(u**2 + 4)
    mu = numerator / denominator
    return mu, u

# Plotting loop over B values
for i, x in enumerate(bx):
    # Create an empty list to store u values for the current B
    # Calculate thetaE for the current mass
   # thetaE_val = thetaE(mass, G, Dsl, Ds, Dl)
    thetaE_val = 1
    mu_array = []
    bxt = []
    by= 0.2 *  thetaE_val
    
    mu, u = compute_mu(bx[i], by, thetaE_val )
    mu_array.append(mu)
    bxt.append(bx[i] / thetaE_val)
    
    
    t = time(Dl, thetaE_val,v)
    td = (mass,Dl,Dsl,Ds,v)
    mu_array  = np.array(mu_array)
    bxt = np.array (bxt)
    
    plt.plot (bxt,  mu_array , 'o')
     


plt.xlabel('Bx/ThetaE')
plt.ylabel('Star Image Magnification factor (u)')
plt.title('Magnification  vs. Impact Parameter x component,  by= 0.2 *  thetaE_val')
plt.show()
sys.exit()
