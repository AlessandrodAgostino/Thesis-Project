from section import section

section(rotation = False,
        N_points = int(snakemake.wildcards.N_points),
        seed = int(snakemake.wildcards.seed),
        saving_path = '',
        n_slices = [-1,0,1])
