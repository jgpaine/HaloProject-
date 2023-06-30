import numpy as np
import matplotlib.pyplot as plt

# Read data from the output file
#data = np.loadtxt('DL=' + str(Dl) + 'Ds=.'+ str(Ds)+'v='+ str(v) + 'TE-' + str(close) +  '.txt', delimiter=',')

data1 = np.loadtxt('PeakMagoutput/DL=5.0Ds=.10.0v=100TE-0.1.txt',delimiter=',')
#data2=np.loadtxt('DL=10.0Ds=25.0v=100TE-0.1.txt',delimiter=',')
data2=np.loadtxt('PeakMagoutput/DL=50.0Ds=550.0v=100TE-0.1.txt',delimiter=',')
data3= np.loadtxt('PeakMagoutput/DL=700.0Ds=710.0v=100TE-0.1.txt',delimiter=',')
data4= np.loadtxt('PeakMagoutput/DL=3.0Ds=503.0v=100TE-0.1.txt',delimiter=',')
data5 = np.loadtxt('PeakMagoutput/DL=5.0Ds=705.0v=100TE-0.1.txt',delimiter=',')

# Extract massarray and magmax arrays from the data
massarray1 = data1[:, 0]
magmax1 = data1[:, 1]
massarray2 = data2[:, 0]
magmax2 = data2[:, 1]
massarray3 = data3[:, 0]
magmax3 = data3[:, 1]
massarray4 = data4[:, 0]
magmax4 = data4[:, 1]
massarray5 = data5[:, 0]
magmax5 = data5[:, 1]

# Plot the data

plt.semilogx(massarray1, magmax1, '.', color='blue', label='DL = 5.0 kpc, Ds = 10 kpc, v=100 km/s, TE= 0.1')
#plt.semilogx(massarray2, magmax2, '.', color='red', label='DL = 10.0 kpc, Ds = 25 kpc, v=100 km/s, TE= 0.1')
plt.semilogx(massarray2, magmax2, '.', color='red', label='DL = 50.0 kpc, Ds = 550 kpc, v=100 km/s, TE= 0.1')
plt.semilogx(massarray3, magmax3, '.', color='green', label='DL = 700.0 kpc, Ds = 710 kpc, v=100 km/s, TE= 0.1')
plt.semilogx(massarray4, magmax4, '.', color='purple', label='DL = 3.0 kpc, Ds = 503 kpc, v=100 km/s, TE= 0.1')
plt.semilogx(massarray5, magmax5, '.', color='orange', label='DL = 5.0 kpc, Ds = 705 kpc, v=100 km/s, TE= 0.1')

plt.xlabel('Point Mass in Solar Masses')
plt.ylabel('Maximum Star Image Magnification factor (u)')
plt.title('Mass vs Maximum Magnification')

#text_message = 'Minimum impact parameter = ' + str(close) + ' * Einstein Radius\n'
#text_message += 'Distance to lens = ' + str(Dl) + ' kpc\n'
#text_message += 'Distance to source = ' + str(Ds) + ' kpc'

#plt.text(1e-6, 2, text_message, fontsize=8, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.legend()
plt.show()

