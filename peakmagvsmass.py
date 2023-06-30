import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# june 15 2023, Make the plot of peak maginfication vs mass of point


#mass = 1  # solar mass
Dl = 0.1  # kps
Dsl =700.0 # kps
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
  #  print(math.sqrt(bx**2 + by**2))
    numerator = u**2 + 2
    denominator = u * math.sqrt(u**2 + 4)
    mu = numerator / denominator
    return mu, u

# make mass array
massarray  = np.logspace(np.log10(1e-6), np.log10(1e9), num=100)
magmax = np.empty(len(massarray))
#print(massarray)

for j, m in enumerate(massarray):
    # Plotting loop over B values
    for i, x in enumerate(bx):
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
    
    
     t = time(Dl, thetaE_val,v)
     td = (massarray[j],Dl,Dsl,Ds,v)
     mu_array  = np.array(mu_array)
     bxt = np.array (bxt)
    
     # plt.plot (bxt,  mu_array , 'o')
    magmax[j] = np.max(mu_array)
    magmax = np.array(magmax)
    #print( magmax[j])
    #plt.plot(massarray[j], magmax[j], 'o')
    plt.semilogx(massarray[j], magmax[j], '.', color='blue')

#create output file
with open('DL=' + str(Dl) + 'Ds='+ str(Ds)+'v='+ str(v) + 'TE-' + str(close) +  '.txt' , 'w') as file:
    for j in range(len(massarray)):
        file.write(f"{massarray[j]}, {magmax[j]}\n")
   


plt.xlabel('Point Mass in Solar Masses')
plt.ylabel('Maxium Star Image Magnification factor (u)')
plt.title('Mass vs Maxium Magnification')

text_message = 'Minium impact parameter = ' + str(close) + ' * Eistein Radius\n'
text_message += 'Distance to lens =' + str(Dl) + ' kpc\n'
text_message += 'Distance to source =' + str(Ds) + ' kpc'

plt.text(1e-6, 2, text_message, fontsize=8, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.show()
sys.exit()
