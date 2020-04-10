import numpy as np
import matplotlib.pyplot as plt

# average = plt.imread(snakemake.input[0])
# noise = plt.imread(snakemake.input[1])

average = plt.imread('../N_6000_seed_42_average.png')
noise = plt.imread('../N_6000_seed_42_noise.png')

average.shape
noise.shape



f_n = 0.8
output = np.divide(average + f_n*noise, 1+f_n)
plt.imsave(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_average+noise.png',output)
