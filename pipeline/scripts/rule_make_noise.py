from section import _draw_noise

fig = _draw_noise()
fig.savefig(f'N_{snakemake.wildcards.N_points}_seed_{snakemake.wildcards.seed}_noise.png', dpi =100)
