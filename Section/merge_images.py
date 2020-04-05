import os
import numpy as np

import matplotlib.pyplot as plt

from section import section, _draw_noise

section(n_slices = 4, rotation = True, N_points = 25000, seed = 290, saving_path = '', noise_density=50)



files = os.listdir('.')
sl_files = [ f for f in files if 'sl' in f]
images = np.array([plt.imread(f) for f in sl_files])
noise = plt.imread(*[ f for f in files if 'noise' in f])

f_n = 0.8
output = np.divide(images.sum(axis=0) + noise*f_n,  4+f_n)
plt.imsave('average_noisy.png', output)
