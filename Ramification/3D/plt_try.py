
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(100)
y = -x

plt.figure()
plt.plot(x)
plt.figure()
plt.plot(y)

for i in plt.get_fignums():
    plt.figure(i)
    plt.savefig('figure%d.png' % i)
