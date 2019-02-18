from uwimg import *

im = load_image("data/Rainier1.png")
detect_and_draw_corners(im, 2, 50, 3)
save_image(im, "corners")

a = load_image("data/Rainier1.png")
b = load_image("data/Rainier2.png")
m = find_and_draw_matches(a, b, 2, 50, 3)
save_image(m, "matches")

im1 = load_image("data/Rainier1.png")
im2 = load_image("data/Rainier2.png")
pan = panorama_image(im1, im2, thresh=50)
save_image(pan, "easy_panorama")

