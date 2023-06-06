import sys
import math
import numpy as np
import matplotlib.pyplot as plt
# june 6 2023 

# Mass in solar masses
#M = np.logspace(0, 4, num=200)
mass = 1*10**20 # solar mass
# Distances in kpc around 10 kps 1 kiloparsec = 3.086 * 10**16 km
Dl = 0.0006 * 3.086 * 10**16
Dsl = 0.0005 * 3.086 * 10**16
Ds = Dl + Dsl
v = 100 # km/s
c2= 9*10**10 #km/s^2

# Gravitational constant in km^2 Mpc M_sun^-1 s^-2
G = 4.301 * 10**-9
# add converstion

# B values in km
bx = [-9,-8,-7,-6,-5,-4,-3,-2,-1,0.000001,1, 2, 3, 4, 5, 6, 7, 8, 9]
length = len(bx)
start_value = 0.1 #by min, closest approach
increment = 0.1
by = np.arange(start_value, start_value + increment * length, increment)


# define this in mass loop
def thetaE(Mass, Grav, Dstol, Dtos, Dtol): # eq.2.70
 
    thetaE = math.sqrt((4 * Grav * Mass * Dstol) / (c2* Dtos * Dtol * Dstol))
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

print(compute_mu(10,10, 2))
print(compute_mu(30,10, 3))

# Plotting loop over B values
for x in range(len(bx)):
    # Create an empty list to store u values for the current B
    # Calculate thetaE for the current mass
    thetaE_val = thetaE(mass, G, Dsl, Ds, Dl) # this angle is too small one the scale of the galaxy for the plot mu=1, in rads
   # print( thetaE_val)
    #   1 Radians = 206264806719.15 Microarcseconds
    theta_uas= thetaE_val * 206264806719.15
    print("theta = ", theta_uas)
    for y in range(len(by)):
      mu_array = []
      u_array = []
      mu, u = compute_mu(bx[x], by[y], theta_uas )
      mu_array.append(mu)
      u_array.append(u)
      
      t = time(Dl, thetaE_val,v) #later check which units for theta
      td = (mass,Dl,Dsl,Ds,v)
      
      mu_array  = np.array(mu_array)
      u_array = np.array(u_array)
      plt.plot(u_array,  mu_array, 'o' )
      print("mu = ", mu)


plt.xlabel('B_total/ThetaE')
plt.ylabel('Star Image Magnification factor (u)')
plt.title('Mag vs. ratio')
plt.show()
sys.exit()
