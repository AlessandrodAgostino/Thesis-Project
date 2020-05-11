from PIL import Image

im = Image.open(snakemake.input[0])
im = im.resize((128,128))
im.save(snakemake.output[0])
