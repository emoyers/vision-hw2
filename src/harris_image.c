#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    return make_image(1,1,1);
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        return copy_image(im);
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    S.data = calloc(im.w*im.h*3, sizeof(float));

    image f_sobel_x = make_gx_filter();
    image f_sobel_y = make_gy_filter();
    image I_x = convolve_image(im, f_sobel_x, 0);
    image I_y = convolve_image(im, f_sobel_y, 0);

    assert( (I_x.h == I_y.h) && (I_x.w == I_y.w) && (I_x.c == I_y.c));

    float pixel_Ix = 0.0f;
    float pixel_Iy = 0.0f;
    float pixel_S = 0.0f;
    int size_column = I_x.w; 
    int size_row = I_x.h;
    int size_channel = I_x.c;

    image IxIx = make_image(I_x.w, I_x.h, I_x.c);
    IxIx.data = calloc(I_x.w*I_x.h*I_x.c, sizeof(float));

    image IyIy = make_image(I_y.w, I_y.h, I_y.c);
    IxIx.data = calloc(I_y.w*I_y.h*I_y.c, sizeof(float));

    image IxIy = make_image(I_x.w, I_x.h, I_x.c);
    IxIy.data = calloc(I_x.w*I_x.h*I_x.c, sizeof(float));

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                pixel_Ix = get_pixel(I_x , i, j ,k);
                pixel_Iy = get_pixel(I_y , i, j ,k);

                set_pixel(IxIx, i, j ,k, pixel_Ix*pixel_Ix);
                set_pixel(IyIy, i, j ,k, pixel_Iy*pixel_Iy);
                set_pixel(IxIy, i, j ,k, pixel_Ix*pixel_Iy);

            }
        }
    }

    image filter = make_gaussian_filter(sigma);
    image IxIx_blur = convolve_image(IxIx, filter, 1);
    image IyIy_blur = convolve_image(IyIy, filter, 1);
    image IxIy_blur = convolve_image(IxIy, filter, 1);

    assert( (IxIx_blur.c == 1) && (IyIy_blur.c == 1) && (IxIy_blur.c == 1));

    size_column = S.w; 
    size_row = S.h;
    size_channel = S.c;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){

                if(k==0)
                {
                    pixel_S = get_pixel(IxIx_blur , i, j ,k);
                }
                else if(k==1)
                {
                    pixel_S = get_pixel(IyIy_blur , i, j ,k);
                }
                else
                {
                    pixel_S = get_pixel(IxIy_blur , i, j ,k);
                }


                set_pixel(S, i, j ,k, pixel_S);

            }
        }
    }

    free_image(f_sobel_x);
    free_image(f_sobel_y);
    free_image(I_x);
    free_image(I_y);
    free_image(IxIx);
    free_image(IyIy);
    free_image(IxIy);
    free_image(filter);
    free_image(IxIx_blur);
    free_image(IyIy_blur);
    free_image(IxIy_blur);
    // TODO: calculate structure matrix for im.
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    R.data = calloc(S.w*S.h*1, sizeof(float));
    float pixel_IxIx = 0.0f;
    float pixel_IyIy = 0.0f;
    float pixel_IxIy = 0.0f;
    float determinant = 0.0f;
    float trace = 0.0f;
    float cornerness_result = 0.0f;
    int size_column = R.w; 
    int size_row = R.h;
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    for(int j = 0; j < size_row; ++j){
        for(int i = 0; i < size_column; ++i){

            pixel_IxIx = get_pixel(S , i, j ,0);
            pixel_IyIy = get_pixel(S , i, j ,1);
            pixel_IxIy = get_pixel(S , i, j ,2);
            determinant = (pixel_IxIx * pixel_IyIy) -(pixel_IxIy * pixel_IxIy);
            trace = (pixel_IxIx + pixel_IyIy);
            cornerness_result = determinant - (0.06 * trace * trace);

            set_pixel(R, i, j ,0, cornerness_result);

        }
    }

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);

    int size_column = r.w; 
    int size_row = r.h;
    int size_channel = r.c;
    float pixel_reference = 0.0f;
    float pixel_comparison = 0.0f;
    int window_size = (2 * w + 1);
    int h_window_size = (int)(window_size / 2);
    int map_pixel_x = 0;
    int map_pixel_y = 0;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){

                pixel_reference = get_pixel(im , i, j ,k);

                /*for neighbors within w:*/
                for(int y = 0; y < window_size; ++y){
                    for(int x = 0; x < window_size; ++x){

                        map_pixel_x = i + x - h_window_size;
                        map_pixel_y = j + y - h_window_size;

                        pixel_comparison = get_pixel(im , map_pixel_x, map_pixel_y ,k);

                        if (pixel_comparison > pixel_reference){

                            r.data[i+j*size_column+k*size_column*size_row] = -999999;
                            /*exit fors x and y*/
                            y =  window_size;
                            x =  window_size;
                        }

                    }
                }
                /*for neighbors within w:*/

            }
        }
    }
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    //TODO: count number of responses over threshold
    int count = 0;

    int size_column = Rnms.w; 
    int size_row = Rnms.h;
    int size_channel = Rnms.c;
    float pixel = 0.0f;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){

                pixel = get_pixel(Rnms, i, j, k);
                if (pixel > thresh){

                    count++;
                }

            }
        }
    }

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int descriptor_count = 0;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){

                pixel = get_pixel(Rnms, i, j, k);
                if ((pixel > thresh) && (descriptor_count < count)){

                    d[descriptor_count].p.x = i;
                    d[descriptor_count].p.y = j;
                    descriptor_count++;
                }

            }
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
