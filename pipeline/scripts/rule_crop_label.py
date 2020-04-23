import numpy as np
import matplotlib.pyplot as plt

images = plt.imread(snakemake.input[0])
crop = (slice(110,710), slice(110,710))
plt.imsave(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_label.png',images[crop])
