from uwimg import *

im = load_image("data/dogbw.png")
s = structure_matrix(im, 2)
c = cornerness_response(s)
nms = nms_image(c,3)
feature_normalize(c)
feature_normalize(s)
feature_normalize(nms)
save_image(s, "dog_structure")
save_image(c, "dog_cornerness")
save_image(nms, "dog_nms")

im = load_image("data/Rainier1.png")
s = structure_matrix(im, 2)
c = cornerness_response(s)
nms = nms_image(c,3)
feature_normalize(c)
feature_normalize(s)
feature_normalize(nms)
save_image(s, "mountain_structure")
save_image(c, "mountain_cornerness")
save_image(nms, "mountain_nms")

im = load_image("data/Rainier1.png")
detect_and_draw_corners(im, 2, 50, 3)
save_image(im, "corners")

