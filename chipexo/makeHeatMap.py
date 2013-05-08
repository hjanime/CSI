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

blue_red = LinearSegmentedColormap("BlueRed", cdict3)

data = np.loadtxt(sys.argv[1])

upper = float(sys.argv[2])


data[data<0.00001] = None
width = data.shape[1] / 2

print width
data[data>upper] = upper

data[:,width:] = data[:,width:] * -1

#data = data.reshape((data.shape[0]*2,data.shape[1]/2))

figure = plt.figure()
plt.imshow(data, aspect='auto', cmap=blue_red)
#plt.imshow(data[:,:width],aspect='auto', cmap=blue_red)
plt.colorbar()
ticks = np.array( [0, width/2-15, width/2, width, width + width/2 , width*1.5 + 15,  2*width] )
labels = [-int(width/2), -15, 0, 'sep', 0,15, int(width/2)]
plt.xticks(ticks, labels)

plt.xlabel('Distance from motif (bp)')
plt.title(sys.argv[1].split('/')[-1].split('.')[0])

tokens = sys.argv[1].split('.')
tokens[-1] = 'heatMap.png'
#plt.savefig('.'.join(tokens))
plt.show()
