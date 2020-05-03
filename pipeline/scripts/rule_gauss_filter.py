from skimage.filters import gaussian
import matplotlib.pyplot as plt
import numpy as np

image = plt.imread(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_average.png')
blur_image = gaussian(image, sigma = 1.5, multichannel = True)
plt.imsave(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_average_blurr.png', blur_image)
