import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# june 15 2023, Make the plot of peak maginfication vs mass of point


#mass = 1  # solar mass
Dl = 5.0  # kps
Dsl = 5.0 # kps
Ds = Dl + Dsl
v = 100 # km/s
c2= 9*10**10 #km/s^2 
# Gravitational constant in km^2 kpc M_sun^-1 s^-2
G = 4.302 * 10**-6
close = 0.1

num_points = 100
num_points_per_side = num_points // 2
log_points = np.logspace(-4, 0, num=num_points_per_side)
#bx = np.concatenate((-log_points[::-1], log_points))
bx = np.concatenate((-log_points[::-1] / 100, log_points / 100))
# bx = [-9,-8,-7,-6,-5,-4,-3,-2,-1,-0.5,-0.25,-0.05,0.000001,0.05,0.25,0.5,1, 2, 3, 4, 5, 6, 7, 8, 9]

def thetaE(Mass, Grav, Dstol, Dtos, Dtol): # eq.2.70
 
    thetaE = math.sqrt((4 * Grav * Mass * Dstol) / (c2 * Dtos * Dtol )) # add c^2
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

# make mass array
massarray  = np.logspace(np.log10(1e-6), np.log10(1e9), num=100)
magmax = np.empty(len(massarray))


for j, m in enumerate(massarray):
    for i, x in enumerate(bx):
    # Plotting loop over B values
     # Calculate thetaE for the current mass
     thetaE_val = thetaE(massarray[j], G, Dsl, Ds, Dl) #radians
     # convert to  arc seconds
     thetaE_as = thetaE_val * 206264
     #thetaE_val = 0.001
     #  1 r =  206264 Arcsecond
     mu_array = []
     bxt = []
     by= close *  thetaE_as
    
     mu, u = compute_mu(bx[i], by, thetaE_as )
     mu_array.append(mu)
     bxt.append(bx[i] / thetaE_as) # leave as _val not _as so bx in rads?
     mu_array  = np.array(mu_array)
     bxt = np.array (bxt)
     
     td = time(Dl,thetaE_as,v)
   # td = duration(massarray[j],Dl,Dsl,Ds,v)
     timeday = []
     timeday = np.array(td)
     
    magmax[j] = np.max(mu_array)
    magmax = np.array(magmax)
    print(timeday)
    timesec = timeday / 86400
    plt.semilogx(massarray[j], timeday, '*')
   
 

plt.xlabel('Point Mass in Solar Masses')
plt.ylabel('Duration in Days')
plt.title('Mass vs Duration')

text_message = 'Minium impact parameter = ' + str(close) + ' * Eistein Radius\n'
text_message += 'Distance to lens =' + str(Dl) + ' kpc\n'
text_message += 'Distance to sorce =' + str(Ds) + ' kpc\n'
text_message += 'Lens Velocity = ' + str(v) + ' km/s^2'

plt.text(1e-5, 0.5, text_message, fontsize=8, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.show()
sys.exit()
