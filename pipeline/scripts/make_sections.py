from section import section

section(rotation = True,
        N_points = int(snakemake.wildcards.N_points),
        seed = int(snakemake.wildcards.seed),
        saving_path = '',
        n_slices = snakemake.config['SLICES'])
