from PIL import Image
import numpy as np
import skfmm

def make_sdf(infile, outfile):
    img = Image.open(infile)
    img.load()
    
    phi = np.asarray(img, dtype="int32")[:,:,0] / 255

    sdf_out = skfmm.distance(-phi, dx = 1 / phi)
    sdf_in = skfmm.distance(phi, dx = 1 / phi)
    sdf_final = (1 + sdf_out - sdf_in) * 0.5
            
    im_sdf = Image.fromarray(sdf_final * 255).convert(img.mode)
    im_sdf.show()
    #im_sdf.save(outfile)

if __name__ == "__main__":
    make_sdf("test_obj/position.png", "test_obj/sdf.png")