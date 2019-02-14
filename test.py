from uwimg import *

im = load_image("data/Rainier1.png")
detect_and_draw_corners(im, 2, 50, 3)
save_image(im, "corners")

a = load_image("data/Rainier1.png")
b = load_image("data/Rainier2.png")
m = find_and_draw_matches(a, b, 2, 50, 3)
save_image(m, "matches")

