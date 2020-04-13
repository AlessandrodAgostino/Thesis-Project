import numpy as np
import matplotlib.pyplot as plt

average = plt.imread(snakemake.input[0])
noise = plt.imread(snakemake.input[1])

crop = (slice(110,710), slice(110,710))

f_n = 0.8
output = np.divide(average[crop] + f_n*noise[crop], 1+f_n)
plt.imsave(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_average+noise.png',output)
