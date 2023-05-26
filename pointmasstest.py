import math
import numpy as np
import matplotlib.pyplot as plt

# Mass in solar masses
M = np.logspace(0, 5, num=100)

# Distances in kpc
Dl = 20.0
Dsl = 10.0
Ds = 30.0
v = 100

# Gravitational constant in km^2 Mpc M_sun^-1 s^-2
G = 4.301 * 10**-9

# B values in km
#B_values = [1, 2, 3, 4, 5, 6, 7, 8, 9]
B_values = [1,1.00001,1.0000002]
def thetaE(Mass, Grav, Dstol, Dtos, Dtol): # eq.2.70
 
    thetaE = math.sqrt((4 * Grav * Mass * Dstol) / (Dtos * Dtol * Dstol))
    return thetaE

def total_magnification(B, thetaE): # eq 4.16
    u = B**2 + 2 * thetaE**2 / (B * np.sqrt(B**2 + 4 * thetaE))
    return u

def time(D,th,velocity): # 5.4 
    t =( (D*th) / velocity ) 
    return t 

# Create an empty list to store u values for different B values
u_values = []

# Plotting loop over B values
for B in B_values:
    # Create an empty list to store u values for the current B
    u_values_b = []
    
    # Loop over mass values
    for mass in M:
        # Calculate thetaE for the current mass
        thetaE_val = thetaE(mass, G, Dsl, Ds, Dl)
        
        # Calculate u for the current B and thetaE
        u = total_magnification(B, thetaE_val)
        
        # Append the u value to the current B's u values list
        u_values_b.append(u)
        
        t = time(Dl, thetaE_val,v)
        print(t)
    
    # Convert the u values list for the current B to a NumPy array
    u_values_b = np.array(u_values_b)
    
    # Plot the mass vs. u for the current B with a different color
    plt.plot(M, u_values_b, 'o', label=f'B = {B}')


    
# Set the plot labels and title
plt.xlabel('Mass (M☉)')
plt.ylabel('Magnification factor (u)')
plt.title('Mass vs. u for Different B Values')
#plt.xscale('log')

# Add a legend to distinguish the different B values
plt.legend()

# Display the plot
plt.show()

