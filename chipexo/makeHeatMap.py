import os,sys
import numpy as np
import pylab as plt
from matplotlib.colors import LinearSegmentedColormap


cdict3 = {'blue':  ((0.0, 0.0, 0.0),
                    (0.25,0.0, 0.0),
                    (0.5, 1.0, 1.0),
                        (0.75,1.0, 1.0),
                    (1.0, 0.4, 1.0)),

    'green': ((0.0, 0.0, 0.0),
        (0.25,0.0, 0.0),
        (0.5, 1.0, 1.0),
        (0.75,0.0, 0.0),
        (1.0, 0.0, 0.0)),

    'red':  ((0.0, 0.0, 0.4),
        (0.25,1.0, 1.0),
        (0.5, 1.0, 1.0),
        (0.75,0.0, 0.0),
        (1.0, 0.0, 0.0)),
    'alpha':  ((0.0, 1.0, 1.0),
               (0.5, 0.3, 0.3),
               (1.0,1.0,1.0))
    }

tokens = sys.argv[1].split('.')
blue_red = LinearSegmentedColormap("BlueRed", cdict3)

data = np.loadtxt(sys.argv[1])

upper = float(sys.argv[2])

width = data.shape[1] / 2

averages = np.average(data, axis=0)
w = np.ones(10)

faves = np.convolve(w/w.sum(), averages[:width], mode='valid')
raves = np.convolve(w/w.sum(), averages[width:], mode='valid')
#averages[width:] = averages[width:] * -1

figure = plt.figure()
#plt.plot(averages[:width])
#plt.plot(averages[width:])
plt.plot(faves)
plt.plot(raves)
tokens[-1] = 'averaged.png'
ticks_half = np.array([0,width/2, width])
labels_half = np.array([-(width/2), 0, width/2])
plt.xticks(ticks_half, labels_half)
plt.savefig('.'.join(tokens), dpi=600)

plt.clf()

#Plot the medians
medians = np.zeros(data.shape[1])
w = np.ones(10)
for i in range(data.shape[1]):
    medians[i] = np.median(data[:,i][data[:,i]>0])

#fmedians = np.convolve(w/w.sum(), medians[:width], mode='valid')
#rmedians = np.convolve(w/w.sum(), medians[width:], mode='valid')
fmedians = medians[:width]
rmedians = medians[width:]

figure = plt.figure()
plt.plot(fmedians)
plt.plot(rmedians)
tokens[-1] = 'median.png'
plt.xticks(ticks_half, labels_half)
plt.savefig('.'.join(tokens), dpi=600)

plt.clf()


data[data<0.00001] = None

print width
data[data>upper] = upper

data[:,width:] = data[:,width:] * -1

#data = data.reshape((data.shape[0]*2,data.shape[1]/2))

figure = plt.figure()
plt.imshow(data, aspect='auto', cmap=blue_red)
#plt.imshow(data[:,:width],aspect='auto', cmap=blue_red)
plt.colorbar()
ticks = np.array( [0,  width/2,  width,  width + width/2 ,  2*width] )
labels = [-int(width/2),  0,  'sep', 0, int(width/2)]
plt.xticks(ticks, labels)

plt.xlabel('Distance from motif (bp)')
plt.title(sys.argv[1].split('/')[-1].split('.')[0])

tokens[-1] = 'heatMap.png'
plt.savefig('.'.join(tokens), dpi=600)
#plt.show()
