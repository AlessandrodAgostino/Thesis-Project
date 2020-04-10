import numpy as np
import matplotlib.pyplot as plt

images = np.array([plt.imread(f) for f in snakemake.input])
output = np.divide(images.sum(axis=0),len(snakemake.input))
plt.imsave(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_average.png',output)
