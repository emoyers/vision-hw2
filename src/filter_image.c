#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    int size_column = im.w; 
    int size_row = im.h;
    int size_channel = im.c;
    float divisor = 0.0f;
    float pixel_value = 0.0f;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                divisor += get_pixel(im, i, j, k);
            }
        }
    }

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                pixel_value = get_pixel(im, i, j, k);
                set_pixel(im, i, j, k, (pixel_value / divisor));
            }
        }
    }

}

image make_box_filter(int w)
{
    // TODO
    image filter_r = make_image(w, w, 1);
    filter_r.data = calloc(w*w*1, sizeof(float));
    float pixel_value = 0.0f;
    int size_column = w; 
    int size_row = w;

    pixel_value = (1.0f / (float)(w * w));

    for(int j = 0; j < size_row; ++j){
        for(int i = 0; i < size_column; ++i){
            set_pixel(filter_r, i, j, 0, pixel_value);
        }
    }

    return filter_r;

}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    int size_column = im.w; 
    int size_row = im.h;
    int size_channel = im.c;
    int size_column_filter = filter.w; 
    int size_row_filter = filter.h;
    int size_channel_filter = filter.c;
    float average_sum = 0.0f;
    int map_pix_filt_2_img_x = 0;
    int map_pix_filt_2_img_y = 0;
    int half_w_filter_x = 0;
    int half_w_filter_y = 0;
    half_w_filter_x = (size_column_filter / 2);
    half_w_filter_y = (size_row_filter / 2);
    assert( (size_channel == size_channel_filter) || (size_channel_filter == 1));
    /*result image*/
    image result_img = make_image(size_column, size_row, size_channel);
    result_img.data = calloc(size_column * size_row * size_channel, sizeof(float));

    /*image loops*/
    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                average_sum = 0.0f;

                /*filter loops*/
                for(int j_f = 0; j_f < size_row_filter; ++j_f){
                    for(int i_f = 0; i_f < size_column_filter; ++i_f){
                        map_pix_filt_2_img_x = i + i_f - half_w_filter_x;
                        map_pix_filt_2_img_y = j + j_f - half_w_filter_y;
                        int k_f = 0;
                        if(size_channel == size_channel_filter){
                            k_f = k;
                        }
                        else{
                            k_f = 0;
                        }
                        average_sum += ((get_pixel(im, map_pix_filt_2_img_x, map_pix_filt_2_img_y, k)) * (get_pixel(filter, i_f, j_f, k_f)));
                    }
                    set_pixel(result_img, i, j, k, average_sum);
                }
                /*filter loops*/
            }
        }
    }
    /*image loops*/

    /*preserve part*/
    size_column = result_img.w; 
    size_row = result_img.h;
    size_channel = result_img.c;
    if(preserve != 1){
        image result_img_n_p = make_image(size_column, size_row, 1);
        result_img_n_p.data = calloc(size_column * size_row * 1, sizeof(float));
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                average_sum = 0.0f;
                for(int k = 0; k < size_channel; ++k){
                    average_sum += get_pixel(result_img, i, j, k);
                }
                set_pixel(result_img_n_p, i, j, 0, average_sum);             
            }
        }    
        result_img = copy_image(result_img_n_p);
        free_image(result_img_n_p);
    }
    return result_img;
}

image make_highpass_filter()
{
    // TODO
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, 0.0f);
    set_pixel(filter, 1, 0, 0, (-1.0f));
    set_pixel(filter, 2, 0, 0, 0.0f);
    set_pixel(filter, 0, 1, 0, (-1.0f));
    set_pixel(filter, 1, 1, 0, (4.0f));
    set_pixel(filter, 2, 1, 0, (-1.0f));
    set_pixel(filter, 0, 2, 0, 0.0f);
    set_pixel(filter, 1, 2, 0, (-1.0f));
    set_pixel(filter, 2, 2, 0, 0.0f);

    return filter;
}

image make_sharpen_filter()
{
    // TODO
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, 0.0f);
    set_pixel(filter, 1, 0, 0, (-1.0f));
    set_pixel(filter, 2, 0, 0, 0.0f);
    set_pixel(filter, 0, 1, 0, (-1.0f));
    set_pixel(filter, 1, 1, 0, (5.0f));
    set_pixel(filter, 2, 1, 0, (-1.0f));
    set_pixel(filter, 0, 2, 0, 0.0f);
    set_pixel(filter, 1, 2, 0, (-1.0f));
    set_pixel(filter, 2, 2, 0, 0.0f);

    return filter;
}

image make_emboss_filter()
{
    // TODO
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, (-2.0f));
    set_pixel(filter, 1, 0, 0, (-1.0f));
    set_pixel(filter, 2, 0, 0, 0.0f);
    set_pixel(filter, 0, 1, 0, (-1.0f));
    set_pixel(filter, 1, 1, 0, 1.0f);
    set_pixel(filter, 2, 1, 0, 1.0f);
    set_pixel(filter, 0, 2, 0, 0.0f);
    set_pixel(filter, 1, 2, 0, 1.0f);
    set_pixel(filter, 2, 2, 0, 2.0f);

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    int w_sigma = (int)(6.0f*sigma);
    if ((w_sigma%2) == 0 ) {
        w_sigma++;
    }
    image gaussian = make_image(w_sigma, w_sigma, 1);
    gaussian.data = calloc(gaussian.w*gaussian.h*1, sizeof(float));
    float gaussian_r = 0.0f;
    int size_column = gaussian.w; 
    int size_row = gaussian.h;
    float x_gaussian = 0.0f;
    float y_gaussian = 0.0f;
    for(int j = 0; j < size_row; ++j){
        for(int i = 0; i < size_column; ++i){
            x_gaussian = (float)(i - (size_column / 2));
            y_gaussian = (float)(j - (size_row / 2));
            gaussian_r = (1.0f/(TWOPI*sigma*sigma))*exp(((-1.0f)*(x_gaussian*x_gaussian+y_gaussian*y_gaussian))/(2*sigma*sigma));
            set_pixel(gaussian, i, j, 0, gaussian_r);             
        }
    }

    l1_normalize(gaussian);
    return gaussian;
}

image add_image(image a, image b)
{
    // TODO
    assert((a.h == b.h) && (a.w == b.w) && (a.c == b.c));
    int size_column = a.w; 
    int size_row = a.h;
    int size_channel = a.c;
    image add_img = make_image(size_column, size_row, size_channel);
    add_img.data = calloc(size_column * size_row * size_channel, sizeof(float));

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                set_pixel(add_img , i , j , k, ( get_pixel(a , i, j ,k) + get_pixel(b , i, j ,k) ) );
            }
        }
    }
    return add_img;
}

image sub_image(image a, image b)
{
    // TODO
    assert((a.h == b.h) && (a.w == b.w) && (a.c == b.c));
    int size_column = a.w; 
    int size_row = a.h;
    int size_channel = a.c;
    image sub_img = make_image(size_column, size_row, size_channel);
    sub_img.data = calloc(size_column * size_row * size_channel, sizeof(float));

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                set_pixel(sub_img , i , j , k, ( get_pixel(a , i, j ,k) - get_pixel(b , i, j ,k) ) );
            }
        }
    }
    return sub_img;
}

image make_gx_filter()
{
    // TODO
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, (-1.0f));
    set_pixel(filter, 1, 0, 0, 0.0f);
    set_pixel(filter, 2, 0, 0, 1.0f);
    set_pixel(filter, 0, 1, 0, (-2.0f));
    set_pixel(filter, 1, 1, 0, 0.0f);
    set_pixel(filter, 2, 1, 0, 2.0f);
    set_pixel(filter, 0, 2, 0, (-1.0f));
    set_pixel(filter, 1, 2, 0, 0.0f);
    set_pixel(filter, 2, 2, 0, 1.0f);

    return filter;
}

image make_gy_filter()
{
    // TODO
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, (-1.0f));
    set_pixel(filter, 1, 0, 0, (-2.0f));
    set_pixel(filter, 2, 0, 0, (-1.0f));
    set_pixel(filter, 0, 1, 0, 0.0f);
    set_pixel(filter, 1, 1, 0, 0.0f);
    set_pixel(filter, 2, 1, 0, 0.0f);
    set_pixel(filter, 0, 2, 0, 1.0f);
    set_pixel(filter, 1, 2, 0, 2.0f);
    set_pixel(filter, 2, 2, 0, 1.0f);

    return filter;
}

void feature_normalize(image im)
{
    // TODO
    int size_column = im.w; 
    int size_row = im.h;
    int size_channel = im.c;
    float min = 0.0f;
    float max = 0.0f;
    float range = 0.0f;
    float pixel_value = 0.0f;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                pixel_value = get_pixel(im , i, j ,k);
                if (pixel_value > max){
                    max = pixel_value;
                }

                if (pixel_value < min)
                {
                    min = pixel_value;
                }
            }
        }
    }

    range = max - min;

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                if (range == 0.0f)
                {
                    set_pixel(im, i, j, k, 0.0f);
                }
                else{
                    pixel_value = get_pixel(im , i, j ,k);
                    pixel_value = ((pixel_value - min) / range);
                    set_pixel(im, i, j, k, pixel_value);
                }

            }
        }
    }



}

image *sobel_image(image im)
{
    // TODO
    image f_sobel_x = make_gx_filter();
    image f_sobel_y = make_gy_filter();
    image im_sobel_x = convolve_image(im, f_sobel_x, 0);
    image im_sobel_y = convolve_image(im, f_sobel_y, 0);
    assert((im_sobel_x.w == im_sobel_y.w) && (im_sobel_x.h == im_sobel_y.h) && (im_sobel_x.c == im_sobel_y.c));

    int size_column = im_sobel_x.w; 
    int size_row = im_sobel_x.h;
    int size_channel = im_sobel_x.c;
    float pixel_sobel_x = 0.0f;
    float pixel_sobel_y = 0.0f;
    float magnitude = 0.0f;
    float angle = 0.0f;

    image *sobel_result;
    sobel_result = calloc(2, sizeof(im_sobel_x));

    image magnitude_sobel = make_image(size_column, size_row, size_channel);
    magnitude_sobel.data = calloc(size_column * size_row * size_channel, sizeof(float));

    image angel_sobel = make_image(size_column, size_row, size_channel);
    angel_sobel.data = calloc(size_column * size_row * size_channel, sizeof(float));

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                pixel_sobel_x = get_pixel(im_sobel_x , i, j ,k);
                pixel_sobel_y = get_pixel(im_sobel_y , i, j ,k);
                magnitude = sqrt((pixel_sobel_x * pixel_sobel_x) + (pixel_sobel_y * pixel_sobel_y));
                set_pixel(magnitude_sobel, i, j, k, magnitude);
                angle = atan2(pixel_sobel_y , pixel_sobel_x);
                set_pixel(angel_sobel, i, j, k, angle);

            }
        }
    }

    sobel_result[0] = magnitude_sobel;
    sobel_result[1] = angel_sobel;

    free_image(f_sobel_x);
    free_image(f_sobel_y);
    free_image(im_sobel_x);
    free_image(im_sobel_y);

    return sobel_result;
}

image colorize_sobel(image im)
{
    // TODO
    image filter = make_gaussian_filter(2);
    image blur_im = convolve_image(im, filter, 1);

    int size_column = blur_im.w; 
    int size_row = blur_im.h;
    int size_channel = blur_im.c;

    image im_colorize_sobel = make_image(size_column, size_row, size_channel);
    im_colorize_sobel.data = calloc(size_column * size_row * size_channel, sizeof(float));

    image *sobel = sobel_image(blur_im);
    feature_normalize(sobel[0]);
    feature_normalize(sobel[1]);

    for(int k = 0; k < size_channel; ++k){
        for(int j = 0; j < size_row; ++j){
            for(int i = 0; i < size_column; ++i){
                if (k==0)
                {
                    set_pixel(im_colorize_sobel, i, j, k, get_pixel(sobel[1], i, j, k));
                }
                else{
                    set_pixel(im_colorize_sobel, i, j, k, get_pixel(sobel[0], i, j, k));
                }

            }
        }
    }

    hsv_to_rgb(im_colorize_sobel);
    free_image (blur_im);
    free_image (filter);
    free(sobel);

    return im_colorize_sobel;
}
