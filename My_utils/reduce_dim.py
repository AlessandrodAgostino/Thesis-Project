from glob import glob
from PIL import Image

"""
Easy peasy tool to reduce images dimensions
"""


for fn in glob('/home/alessandro/Pictures/super_ris/reduced/*.jpg'):
    im = Image.open(fn)
    im = im.resize((25,25))
    # dir_pos = '/'.join(fn.split('/')[:-2])
    # fn_no_ex = im_name = fn.split('/')[-1].split('.')[0]
    # im.save(dir_pos + '/reduced' + '/' + fn_no_ex + '.jpg')
    im.save(fn)
