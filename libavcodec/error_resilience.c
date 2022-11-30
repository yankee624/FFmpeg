/*
 * Error resilience / concealment
 *
 * Copyright (c) 2002-2004 Michael Niedermayer <michaelni@gmx.at>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * Error resilience / concealment.
 */

#include <limits.h>

#include "libavutil/internal.h"
#include "avcodec.h"
#include "error_resilience.h"
#include "me_cmp.h"
#include "mpegutils.h"
#include "mpegvideo.h"
#include "rectangle.h"
#include "threadframe.h"

/**
 * @param stride the number of MVs to get to the next row
 * @param mv_step the number of MVs per row or column in a macroblock
 */
static void set_mv_strides(ERContext *s, ptrdiff_t *mv_step, ptrdiff_t *stride)
{
    if (s->avctx->codec_id == AV_CODEC_ID_H264) {
        av_assert0(s->quarter_sample);
        *mv_step = 4;
        *stride  = s->mb_width * 4;
    } else {
        *mv_step = 2;
        *stride  = s->b8_stride;
    }
}

/**
 * Replace the current MB with a flat dc-only version.
 */
static void put_dc(ERContext *s, uint8_t *dest_y, uint8_t *dest_cb,
                   uint8_t *dest_cr, int mb_x, int mb_y)
{
    int *linesize = s->cur_pic.f->linesize;
    int dc, dcu, dcv, y, i;
    for (i = 0; i < 4; i++) {
        dc = s->dc_val[0][mb_x * 2 + (i &  1) + (mb_y * 2 + (i >> 1)) * s->b8_stride];
        if (dc < 0)
            dc = 0;
        else if (dc > 2040)
            dc = 2040;
        for (y = 0; y < 8; y++) {
            int x;
            for (x = 0; x < 8; x++)
                dest_y[x + (i &  1) * 8 + (y + (i >> 1) * 8) * linesize[0]] = dc / 8;
        }
    }
    dcu = s->dc_val[1][mb_x + mb_y * s->mb_stride];
    dcv = s->dc_val[2][mb_x + mb_y * s->mb_stride];
    if (dcu < 0)
        dcu = 0;
    else if (dcu > 2040)
        dcu = 2040;
    if (dcv < 0)
        dcv = 0;
    else if (dcv > 2040)
        dcv = 2040;

    if (dest_cr)
    for (y = 0; y < 8; y++) {
        int x;
        for (x = 0; x < 8; x++) {
            dest_cb[x + y * linesize[1]] = dcu / 8;
            dest_cr[x + y * linesize[2]] = dcv / 8;
        }
    }
}

static void filter181(int16_t *data, int width, int height, ptrdiff_t stride)
{
    int x, y;

    /* horizontal filter */
    for (y = 1; y < height - 1; y++) {
        int prev_dc = data[0 + y * stride];

        for (x = 1; x < width - 1; x++) {
            int dc;
            dc = -prev_dc +
                 data[x     + y * stride] * 8 -
                 data[x + 1 + y * stride];
            dc = (av_clip(dc, INT_MIN/10923, INT_MAX/10923 - 32768) * 10923 + 32768) >> 16;
            prev_dc = data[x + y * stride];
            data[x + y * stride] = dc;
        }
    }

    /* vertical filter */
    for (x = 1; x < width - 1; x++) {
        int prev_dc = data[x];

        for (y = 1; y < height - 1; y++) {
            int dc;

            dc = -prev_dc +
                 data[x +  y      * stride] * 8 -
                 data[x + (y + 1) * stride];
            dc = (av_clip(dc, INT_MIN/10923, INT_MAX/10923 - 32768) * 10923 + 32768) >> 16;
            prev_dc = data[x + y * stride];
            data[x + y * stride] = dc;
        }
    }
}
/**
 * Replace the current MB with a flat dc-only version.
 */
static void put_dc_edge(ERContext *s, uint8_t *dest_y, uint8_t *dest_cb,
                   uint8_t *dest_cr, int mb_x, int mb_y, int mb_stride)
{
    int *linesize = s->cur_pic.f->linesize;

    int edge_idx = 0;
    int mb_xy = mb_y * s->mb_stride + mb_x;
    int x, y;
    float x_inc = 0;
    float y_inc = 0;
    int reference_pixels[4] = {0,0,0,0};

    float curr_x, curr_y;
    int curr_pix_x, curr_pix_y;
    int curr_mb_x , curr_mb_y, curr_mb_xy;
    int error_j = 0;
    float dist_right, dist_left;
    int right_x, right_y, left_x, left_y;


   

    float direction = tan(s->cur_pic.mb_edge_dir[mb_y * mb_stride + mb_x]);
    printf("dir: %f\n", s->cur_pic.mb_edge_dir[mb_y * mb_stride + mb_x]*180/3.14);
    if(abs(direction) >= 1){
        x_inc = 1 / direction;
        y_inc = 1;
    } else{
        x_inc = 1;
        y_inc = direction;
    }


    for (y = 0; y < 16; y++) {
        for (x = 0; x < 16; x++){

            curr_x = (float) (x + mb_x * 16);
            curr_y = (float) (y + mb_y * 16);
            
            while(1){
                curr_x += x_inc;
                curr_y -= y_inc;

                curr_pix_x = round(curr_x);
                curr_pix_y = round(curr_y);

                curr_mb_x = curr_pix_x / 16;
                curr_mb_y = curr_pix_y / 16;

                curr_mb_xy = curr_mb_x + curr_mb_y * s->mb_stride;
                int error_j= s->error_status_table[curr_mb_xy];
                if(!(error_j&ER_DC_ERROR)){
                    reference_pixels[0] = curr_pix_x;
                    reference_pixels[1] = curr_pix_y;
                    break;
                }

                if(curr_pix_x >= s->cur_pic.f->width || curr_pix_x < 0 || curr_pix_y >= s->cur_pic.f->height || curr_pix_y < 0) break;
            }
            while(1){
                curr_x -= x_inc;
                curr_y += y_inc;

                curr_pix_x = round(curr_x);
                curr_pix_y = round(curr_y);

                curr_mb_x = curr_pix_x / 16;
                curr_mb_y = curr_pix_y / 16;

                curr_mb_xy = curr_mb_x + curr_mb_y * s->mb_stride;
                int error_j= s->error_status_table[curr_mb_xy];
                if(!(error_j&ER_DC_ERROR)){
                    reference_pixels[2] = curr_pix_x;
                    reference_pixels[3] = curr_pix_y;
                    break;
                }

                if(curr_pix_x >= s->cur_pic.f->width || curr_pix_x < 0 || curr_pix_y >= s->cur_pic.f->height || curr_pix_y < 0) break;
            }

            right_x = reference_pixels[0];
            right_y = reference_pixels[1];
            left_x = reference_pixels[2];
            left_y = reference_pixels[3];
        
            dist_right = sqrt((x-right_x)*(x-right_x)+(y-right_y)*(y-right_y));
            dist_left = sqrt((x-left_x)*(x-left_x) + (y-left_y)*(y-left_y));
          
            if(x + mb_x * 16 > s->cur_pic.f->width || y + mb_y * 16 > s->cur_pic.f->height) continue;

            //printf("%d %d %d %d %d %d %f\n",right_x, right_y, left_x, left_y, dest_y[left_x + left_y * linesize[0]], dest_y[right_x + right_y * linesize[0]], s->cur_pic.mb_edge_dir[mb_y * mb_stride + mb_x]*180/3.14);
            //dest_y[x + mb_x * 16 + (y + mb_y * 16) * linesize[0]] = round((dist_right * dest_y[left_x + left_y * linesize[0]] + dist_left * dest_y[right_x + right_y * linesize[0]]) / (dist_right + dist_left));
           // dest_y[x + y * linesize[0]] = round((dist_right * dest_y[left_x + left_y * linesize[0]] + dist_left * dest_y[right_x + right_y * linesize[0]]) / (dist_right + dist_left));
            dest_y[x + mb_x * 16 + (y + mb_y * 16) * linesize[0]] = round((dist_right * dest_y[left_x + left_y * linesize[0]] + dist_left * dest_y[right_x + right_y * linesize[0]]) / (dist_right + dist_left));
            dest_cb[(y + mb_y * 16)/2*linesize[1] + (x + mb_x * 16)/2] = round((dist_right * dest_cb[left_y/2*linesize[1] + left_x/2] + dist_left * dest_cb[right_y/2*linesize[1] + right_x/2]) / (dist_right + dist_left));
            dest_cr[(y + mb_y * 16)/2*linesize[2] + (x + mb_x * 16)/2] = round((dist_right * dest_cr[left_y/2*linesize[2] + left_x/2] + dist_left * dest_cr[right_y/2*linesize[2] + right_x/2]) / (dist_right + dist_left));
       }
            
    
    }



    /*for (y = 0; y < 8; y++) {
        for (x = 0; x < 8; x++){
            if(x < 4 && y < 4 && x >= y) edge_idx = 0;
            else if(x >= 4 && y < 4 && (x-4) < 4-y) edge_idx = 1;
            else if(x < 4 && y >= 4 && x >= 7-y) edge_idx = 2;
            else if(x >= 4 && y >= 4 && x <= y) edge_idx = 3;
            else if(x < 4 && y < 4 && x >= y) edge_idx = 4;
            else if(x < 4 && y >= 4 && x < 7-y) edge_idx = 5;
            else if(x >= 4 && y < 4 && (x-4) >= 4-y) edge_idx = 6;
            else if(x >= 4 && y >= 4 && x > y) edge_idx = 7;
        }
            
    
    }*/

}

/**
 * guess the dominant edges of all macroblocks
 * dominant edges are 8 values two values per one side of each macroblock
 * right-top/right-bottom/left-top/left-bottom/bottom-left/bottom-right/top-left/top-right
 * this version uses edge information 
 * @param w     width in 8 pixel blocks
 * @param h     height in 8 pixel blocks
*/
static void guess_dominant_edge(ERContext *s, int w, int h, ptrdiff_t stride)
{
    int b_x, b_y;

    float mag[8];
    float dir[8];

    // left to right
    for(b_y=0; b_y<h; b_y++){
        for(b_x=0; b_x<w; b_x++){
            int mb_index_j= b_x + b_y*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_AC_ERROR)){
                mag[6] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 6];
                mag[7] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 7];
                dir[6] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 6];
                dir[7] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 7];
                mag[4] = mag[6]; mag[5] = mag[7];
                dir[4] = dir[6]; dir[5] = dir[7];
               
            }
            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 4] < mag[4]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 4] = mag[4];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 4] = dir[4];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 5] < mag[5]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 5] = mag[5];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 5] = dir[5];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 6] < mag[6]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 6] = mag[6];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 6] = dir[6];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 7] < mag[7]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 7] = mag[7];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 7] = dir[7];
            }
            
            
        }
        
        //right to left
        for(b_x=w-1; b_x>=0; b_x--){
            int mb_index_j= b_x + b_y*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_AC_ERROR)){
                mag[4] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 4];
                mag[5] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 5];
                dir[4] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 4];
                dir[5] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 5];
                mag[6] = mag[4]; mag[7] = mag[5];
                dir[6] = dir[4]; dir[7] = dir[5];
               
            }
            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 4] < mag[4]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 4] = mag[4];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 4] = dir[4];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 5] < mag[5]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 5] = mag[5];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 5] = dir[5];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 6] < mag[6]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 6] = mag[6];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 6] = dir[6];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 7] < mag[7]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 7] = mag[7];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 7] = dir[7];
            }
        }
    }
    for(b_x=0; b_x<w; b_x++){
        
        for(b_y=0; b_y<h; b_y++){
            int mb_index_j= b_x + b_y*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_AC_ERROR)){
                mag[2] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 2];
                mag[3] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 3];
                dir[2] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 2];
                dir[3] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 3];
                mag[0] = mag[2]; mag[1] = mag[3];
                dir[0] = dir[2]; dir[1] = dir[3];
               
            }
            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 0] < mag[0]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 0] = mag[0];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 0] = dir[0];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 1] < mag[1]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 1] = mag[1];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 1] = dir[1];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 2] < mag[2]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 2] = mag[2];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 2] = dir[2];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 3] < mag[3]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 3] = mag[3];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 3] = dir[3];
            }
            
        }
        
        for(b_y=h-1; b_y>=0; b_y--){
            int mb_index_j= b_x + b_y*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_AC_ERROR)){
                mag[0] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 0];
                mag[1] = s->cur_pic.mb_edge_mag[mb_index_j * 8 + 1];
                dir[0] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 0];
                dir[1] = s->cur_pic.mb_edge_dir[mb_index_j * 8 + 1];
                mag[2] = mag[0]; mag[3] = mag[1];
                dir[2] = dir[0]; dir[3] = dir[1];
               
            }
            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 0] < mag[0]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 0] = mag[0];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 0] = dir[0];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 1] < mag[1]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 1] = mag[1];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 1] = dir[1];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 2] < mag[2]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 2] = mag[2];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 2] = dir[2];
            }

            if(s->cur_pic.mb_edge_mag[mb_index_j * 8 + 3] < mag[3]) {
                s->cur_pic.mb_edge_mag[mb_index_j * 8 + 3] = mag[3];
                s->cur_pic.mb_edge_dir[mb_index_j * 8 + 3] = dir[3];
            }
        }
    }

}





static int convolution(ERContext *s, int kernel[3][3], int row, int col) {
	int i, j, sum = 0;
    int *linesize = s->cur_pic.f->linesize;
    int mb_stride = s->mb_stride;

    int curr_x, curr_y;
    int curr_mb_x, curr_mb_y, curr_mb_xy;
    int error_j;
    
	for (i = -1; i < 2; i++) {
		for (j = -1; j < 2; j++) {
            curr_x = j + col;
            curr_y = i + row;
            curr_mb_x = curr_x / 16;
            curr_mb_y = curr_y / 16;
            curr_mb_xy = curr_mb_x + curr_mb_y * mb_stride;
            error_j = s->error_status_table[curr_mb_xy];
            
            if(error_j & ER_MB_ERROR){
                curr_x = col;
                curr_y = row - i;
                curr_mb_x = curr_x / 16;
                curr_mb_y = curr_y / 16;
                curr_mb_xy = curr_mb_x + curr_mb_y * mb_stride;
                error_j = s->error_status_table[curr_mb_xy];
                if(!(error_j & ER_MB_ERROR)){
                    sum += s->cur_pic.f->data[0][curr_y*linesize[0]+curr_x] * kernel[i+1][j+1];
                    continue;
                }

                curr_x = col - j;
                curr_y = row;
                curr_mb_x = curr_x / 16;
                curr_mb_y = curr_y / 16;
                curr_mb_xy = curr_mb_x + curr_mb_y * mb_stride;
                error_j = s->error_status_table[curr_mb_xy];
                if(!(error_j & ER_MB_ERROR)){
                    sum += s->cur_pic.f->data[0][curr_y*linesize[0]+curr_x] * kernel[i+1][j+1];
                    continue;
                }

            } 
            sum += s->cur_pic.f->data[0][curr_y*linesize[0]+curr_x] * kernel[i+1][j+1];

		}
	}
	return sum;
}

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
static void edge_detection_mb(ERContext *s, int mb_x, int mb_y){
    int width = s->cur_pic.f->width; int height = s->cur_pic.f->height;
    
    int start_x = max(16*mb_x, 1);
    int end_x = min(16*(mb_x+1), width - 1);
    int start_y = max(16*mb_y, 1);
    int end_y = min(16*(mb_y+1), height - 1);

    int my[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};
	int mx[3][3] = {
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1}
	};
    int i, j, gx, gy;
    for (i = start_y; i < end_y; i++) {
		for (j = start_x; j < end_x; j++) {
			gx = convolution(s, mx, i, j);
			gy = convolution(s, my, i, j);
            s->cur_pic.edge_mag[i * width + j] = sqrt(gx*gx + gy*gy);
            s->cur_pic.edge_dir[i * width + j] = atan2(gy, gx);
		}
	}

}

static void calculate_dominant_edge(ERPicture p, int mb_x, int mb_y, int width, int height, float* dst_dir, float* dst_mag)
{
    float max_mag = 0;
    float max_dir = 0;
    float mag = 0;
    float mag_sum = 0;
    float dir_sum = 0;
    int curr_x, curr_y;
    
    

}
#define LEFT 0
#define RIGHT 1
#define TOP 2
#define BOTTOM 3
// get edge direction with maximum magnitude 
// for each lossy mb, get closest macroblock with no error from 4 direction (left, right, bottom, top) 
// from the closest macroblock get the biggest edge direction from the boundary pixels
// example) left macroblock, look at the pixels at the right side boundary, get the largest edge direction
static void get_dominant_edge(ERContext *s, int mb_width_idx, int mb_height_idx, int direction, float *result)
{
    float max_magnitude = 0;
    float max_direction = 0;
    int width = s->cur_pic.f->width;
    if(direction == LEFT){
        for(int i = 0 ; i < 16 ; i++){
            if(s->cur_pic.edge_mag[(mb_height_idx * 16 + i) * width + mb_width_idx * 16] > max_magnitude){
                max_magnitude = s->cur_pic.edge_mag[(mb_height_idx * 16 + i) * width + mb_width_idx * 16];
                max_direction = s->cur_pic.edge_dir[(mb_height_idx * 16 + i) * width + mb_width_idx * 16];
            }
        }
    } else if(direction == RIGHT){
        for(int i = 0 ; i < 16 ; i++){
            if(s->cur_pic.edge_mag[(mb_height_idx * 16 + i) * width + (mb_width_idx + 1) * 16 - 1] > max_magnitude){
                max_magnitude = s->cur_pic.edge_mag[(mb_height_idx * 16 + i) * width + (mb_width_idx + 1) * 16 - 1];
                max_direction = s->cur_pic.edge_dir[(mb_height_idx * 16 + i) * width + (mb_width_idx + 1) * 16 - 1];
            }
        }
    } else if(direction == TOP){
        for(int i = 0 ; i < 16 ; i++){
            if(s->cur_pic.edge_mag[mb_height_idx * 16 * width + mb_width_idx  * 16 + i] > max_magnitude){
                max_magnitude = s->cur_pic.edge_mag[mb_height_idx * 16 * width + mb_width_idx  * 16 + i];
                max_direction = s->cur_pic.edge_dir[mb_height_idx * 16 * width + mb_width_idx  * 16 + i];
            }
            
        }

    } else if(direction == BOTTOM){
        for(int i = 0 ; i < 16 ; i++){
            if(s->cur_pic.edge_mag[((mb_height_idx+1)* 16 - 1) * width + mb_width_idx  * 16 + i] > max_magnitude){
                max_magnitude = s->cur_pic.edge_mag[((mb_height_idx+1)* 16 - 1) * width + mb_width_idx  * 16 + i];
                max_direction = s->cur_pic.edge_dir[((mb_height_idx+1)* 16 - 1) * width + mb_width_idx  * 16 + i];
            }
        }
    }
    result[0] = max_magnitude;
    result[1] = max_direction;
}

static void edge_detection_picture(ERContext *s, int mb_width, int mb_height, int mb_stride)
{
    int width = s->cur_pic.f->width; int height = s->cur_pic.f->height;
    s->cur_pic.edge_mag = av_malloc_array(width, height * sizeof(float));
    s->cur_pic.edge_dir = av_malloc_array(width, height * sizeof(float));
    s->cur_pic.mb_edge_dir = av_malloc_array(mb_stride, mb_height * sizeof(float));
    s->cur_pic.mb_edge_mag = av_malloc_array(mb_stride, mb_height * sizeof(float));
  
    int i, j, gx, gy;
	int my[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};
	int mx[3][3] = {
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1}
	};


    int curr_mb_x, curr_mb_y, curr_mb_xy;
    int error_j;
	
	for (i = 1; i < height - 1; i++) {
		for (j = 1; j < width - 1; j++) {
            curr_mb_x = j / 16;
            curr_mb_y = i / 16;
            curr_mb_xy = curr_mb_x + curr_mb_y * mb_stride;
            error_j = s->error_status_table[curr_mb_xy];
            if(error_j & ER_MB_ERROR){
                s->cur_pic.edge_mag[i * width + j] = 0;
                s->cur_pic.edge_dir[i * width + j] = 0;
                continue;
            }


			gx = convolution(s, mx, i, j);
			gy = convolution(s, my, i, j);
            s->cur_pic.edge_mag[i * width + j] = sqrt(gx*gx + gy*gy);
            s->cur_pic.edge_dir[i * width + j] = atan2(gy, gx);
		}
	}


    float curr_mag, curr_dir;
    float curr_results[2];
    for(int i = 0 ; i < mb_height ; i++){
        curr_mag = 0;
        curr_dir = 0;
        for(int j = 0 ; j < mb_width ; j++){
            error_j = s->error_status_table[i * mb_stride + j];
            if(!(error_j & ER_MB_ERROR ) ){
                get_dominant_edge(s, j, i, RIGHT, curr_results);
                curr_mag = curr_results[0];
                curr_dir = curr_results[1];
            }

            if(curr_mag > s->cur_pic.mb_edge_mag[i*mb_stride + j]){
                s->cur_pic.mb_edge_mag[i*mb_stride + j] = curr_mag;
                s->cur_pic.mb_edge_dir[i*mb_stride + j] = curr_dir;
            }

        }
    }
    for(int i = 0 ; i < mb_height ; i++){
        curr_mag = 0;
        curr_dir = 0;
        for(int j = mb_width-1 ; j >=0 ; j--){

            error_j = s->error_status_table[i * mb_stride + j];
            if(!(error_j & ER_MB_ERROR )) {
                get_dominant_edge(s, j, i, LEFT, curr_results);
                curr_mag = curr_results[0];
                curr_dir = curr_results[1];
            }

            if(curr_mag > s->cur_pic.mb_edge_mag[i*mb_stride + j]){
                s->cur_pic.mb_edge_mag[i*mb_stride + j] = curr_mag;
                s->cur_pic.mb_edge_dir[i*mb_stride + j] = curr_dir;
            }
        }
    }
    for(int j = 0 ; j < mb_width ; j++){
        curr_mag = 0;
        curr_dir = 0;
        for(int i = 0 ; i < mb_height ; i++){
            error_j = s->error_status_table[i * mb_stride + j];
            if(!(error_j & ER_MB_ERROR )) {
                get_dominant_edge(s, j, i, BOTTOM, curr_results);
                curr_mag = curr_results[0];
                curr_dir = curr_results[1];
            }

            if(curr_mag > s->cur_pic.mb_edge_mag[i*mb_stride + j]){
                s->cur_pic.mb_edge_mag[i*mb_stride + j] = curr_mag;
                s->cur_pic.mb_edge_dir[i*mb_stride + j] = curr_dir;
            }

        }
    }
    for(int j = 0 ; j < mb_width ; j++){
        curr_mag = 0;
        curr_dir = 0;
        for(int i = mb_height - 1 ; i >=0 ; i--){
            error_j = s->error_status_table[i * mb_stride + j];
            if(!(error_j & ER_MB_ERROR )) {
                get_dominant_edge(s, j, i, TOP, curr_results);
                curr_mag = curr_results[0];
                curr_dir = curr_results[1];
            }

            if(curr_mag > s->cur_pic.mb_edge_mag[i*mb_stride + j]){
                s->cur_pic.mb_edge_mag[i*mb_stride + j] = curr_mag;
                s->cur_pic.mb_edge_dir[i*mb_stride + j] = curr_dir;
            }

        }
    }

}

/*static void edge_detection(ERContext *s)
{

    edge_detection_picture(s->cur_pic); 
    if(s->last_pic.f != NULL){
        edge_detection_picture(s->last_pic);
    }
    av_free(s->cur_pic.edge_dir);
    av_free(s->cur_pic.edge_mag);
    if(s->last_pic.f != NULL){
        av_free(s->last_pic.edge_dir);
        av_free(s->last_pic.edge_mag);
    }
    s->cur_pic.edge_mag = av_malloc_array(s->cur_pic.f->width, s->cur_pic.f->height*sizeof(float));
    s->cur_pic.edge_dir = av_malloc_array(s->cur_pic.f->width, s->cur_pic.f->height*sizeof(float));
    s->cur_pic.mb_edge_mag = av_malloc_array(s->mb_stride, s->mb_height*sizeof(float));
    s->cur_pic.mb_edge_dir = av_malloc_array(s->mb_stride, s->mb_height*sizeof(float));


    int i, j, gx, gy;
	int my[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};
	int mx[3][3] = {
		{-1, -2, -1},
		{0, 0, 0},
		{1, 2, 1}
	};
	
	for (i = 1; i < s->cur_pic.f->height - 2; i++) {
		for (j = 1; j < s->cur_pic.f->width - 2; j++) {
			gx = convolution(s->cur_pic.f, mx, i, j);
			gy = convolution(s->cur_pic.f, my, i, j);
            s->cur_pic.edge_mag[i * s->cur_pic.f->width + j] = sqrt(gx*gx + gy*gy);
            s->cur_pic.edge_dir[i * s->cur_pic.f->width + j] = atan2(gy, gx);
		}
	}
}*/

/**
 * guess the dc of blocks which do not have an undamaged dc
 * @param w     width in 8 pixel blocks
 * @param h     height in 8 pixel blocks
 */
static void guess_dc(ERContext *s, int16_t *dc, int w,
                     int h, ptrdiff_t stride, int is_luma)
{

    int b_x, b_y;
    int16_t  (*col )[4] = av_malloc_array(stride, h*sizeof( int16_t)*4);
    uint32_t (*dist)[4] = av_malloc_array(stride, h*sizeof(uint32_t)*4);

    if(!col || !dist) {
        av_log(s->avctx, AV_LOG_ERROR, "guess_dc() is out of memory\n");
        goto fail;
    }


    for(b_y=0; b_y<h; b_y++){
        int color= 1024;
        int distance= -1;
        for(b_x=0; b_x<w; b_x++){
            int mb_index_j= (b_x>>is_luma) + (b_y>>is_luma)*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_DC_ERROR)){
                color= dc[b_x + b_y*stride];
                distance= b_x;
            }
            col [b_x + b_y*stride][1]= color;
            dist[b_x + b_y*stride][1]= distance >= 0 ? b_x-distance : 9999;
        }
        color= 1024;
        distance= -1;
        for(b_x=w-1; b_x>=0; b_x--){
            int mb_index_j= (b_x>>is_luma) + (b_y>>is_luma)*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_DC_ERROR)){
                color= dc[b_x + b_y*stride];
                distance= b_x;
            }
            col [b_x + b_y*stride][0]= color;
            dist[b_x + b_y*stride][0]= distance >= 0 ? distance-b_x : 9999;
        }
    }
    for(b_x=0; b_x<w; b_x++){
        int color= 1024;
        int distance= -1;
        for(b_y=0; b_y<h; b_y++){
            int mb_index_j= (b_x>>is_luma) + (b_y>>is_luma)*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_DC_ERROR)){
                color= dc[b_x + b_y*stride];
                distance= b_y;
            }
            col [b_x + b_y*stride][3]= color;
            dist[b_x + b_y*stride][3]= distance >= 0 ? b_y-distance : 9999;
        }
        color= 1024;
        distance= -1;
        for(b_y=h-1; b_y>=0; b_y--){
            int mb_index_j= (b_x>>is_luma) + (b_y>>is_luma)*s->mb_stride;
            int error_j= s->error_status_table[mb_index_j];
            int intra_j = IS_INTRA(s->cur_pic.mb_type[mb_index_j]);
            if(intra_j==0 || !(error_j&ER_DC_ERROR)){
                color= dc[b_x + b_y*stride];
                distance= b_y;
            }
            col [b_x + b_y*stride][2]= color;
            dist[b_x + b_y*stride][2]= distance >= 0 ? distance-b_y : 9999;
        }
    }

    for (b_y = 0; b_y < h; b_y++) {
        for (b_x = 0; b_x < w; b_x++) {
            int mb_index, error, j;
            int64_t guess, weight_sum;
            mb_index = (b_x >> is_luma) + (b_y >> is_luma) * s->mb_stride;
            error    = s->error_status_table[mb_index];

            if (IS_INTER(s->cur_pic.mb_type[mb_index]))
                continue; // inter
            if (!(error & ER_DC_ERROR))
                continue; // dc-ok

            weight_sum = 0;
            guess      = 0;
            for (j = 0; j < 4; j++) {
                int64_t weight  = 256 * 256 * 256 * 16 / FFMAX(dist[b_x + b_y*stride][j], 1);
                guess          += weight*(int64_t)col[b_x + b_y*stride][j];
                weight_sum     += weight;
            }
            guess = (guess + weight_sum / 2) / weight_sum;
            dc[b_x + b_y * stride] = guess;
        }
    }

fail:
    av_freep(&col);
    av_freep(&dist);
}

/**
 * simple horizontal deblocking filter used for error resilience
 * @param w     width in 8 pixel blocks
 * @param h     height in 8 pixel blocks
 */

 ////  h_block_filter(s, s->cur_pic.f->data[0], s->mb_width * 2,
     //                  s->mb_height * 2, linesize[0], 1);
static void h_block_filter(ERContext *s, uint8_t *dst, int w,
                           int h, ptrdiff_t stride, int is_luma)
{
    int b_x, b_y;
    ptrdiff_t mvx_stride, mvy_stride;
    const uint8_t *cm = ff_crop_tab + MAX_NEG_CROP;
    set_mv_strides(s, &mvx_stride, &mvy_stride);
    mvx_stride >>= is_luma;
    mvy_stride *= mvx_stride;

    for (b_y = 0; b_y < h; b_y++) {
        for (b_x = 0; b_x < w - 1; b_x++) {
            int y;
            int left_status  = s->error_status_table[( b_x      >> is_luma) + (b_y >> is_luma) * s->mb_stride];
            int right_status = s->error_status_table[((b_x + 1) >> is_luma) + (b_y >> is_luma) * s->mb_stride];
            int left_intra   = IS_INTRA(s->cur_pic.mb_type[( b_x      >> is_luma) + (b_y >> is_luma) * s->mb_stride]);
            int right_intra  = IS_INTRA(s->cur_pic.mb_type[((b_x + 1) >> is_luma) + (b_y >> is_luma) * s->mb_stride]);
            int left_damage  = left_status & ER_MB_ERROR;
            int right_damage = right_status & ER_MB_ERROR;
            int offset       = b_x * 8 + b_y * stride * 8;
            int16_t *left_mv  = s->cur_pic.motion_val[0][mvy_stride * b_y + mvx_stride *  b_x];
            int16_t *right_mv = s->cur_pic.motion_val[0][mvy_stride * b_y + mvx_stride * (b_x + 1)];
            if (!(left_damage || right_damage))
                continue; // both undamaged
            if ((!left_intra) && (!right_intra) &&
                FFABS(left_mv[0] - right_mv[0]) +
                FFABS(left_mv[1] + right_mv[1]) < 2)
                continue;

            for (y = 0; y < 8; y++) {
                int a, b, c, d;

                a = dst[offset + 7 + y * stride] - dst[offset + 6 + y * stride];
                b = dst[offset + 8 + y * stride] - dst[offset + 7 + y * stride];
                c = dst[offset + 9 + y * stride] - dst[offset + 8 + y * stride];

                d = FFABS(b) - ((FFABS(a) + FFABS(c) + 1) >> 1);
                d = FFMAX(d, 0);
                if (b < 0)
                    d = -d;

                if (d == 0)
                    continue;

                if (!(left_damage && right_damage))
                    d = d * 16 / 9;

                if (left_damage) {
                    dst[offset + 7 + y * stride] = cm[dst[offset + 7 + y * stride] + ((d * 7) >> 4)];
                    dst[offset + 6 + y * stride] = cm[dst[offset + 6 + y * stride] + ((d * 5) >> 4)];
                    dst[offset + 5 + y * stride] = cm[dst[offset + 5 + y * stride] + ((d * 3) >> 4)];
                    dst[offset + 4 + y * stride] = cm[dst[offset + 4 + y * stride] + ((d * 1) >> 4)];
                }
                if (right_damage) {
                    dst[offset + 8 + y * stride] = cm[dst[offset +  8 + y * stride] - ((d * 7) >> 4)];
                    dst[offset + 9 + y * stride] = cm[dst[offset +  9 + y * stride] - ((d * 5) >> 4)];
                    dst[offset + 10+ y * stride] = cm[dst[offset + 10 + y * stride] - ((d * 3) >> 4)];
                    dst[offset + 11+ y * stride] = cm[dst[offset + 11 + y * stride] - ((d * 1) >> 4)];
                }
            }
        }
    }
}

/**
 * simple vertical deblocking filter used for error resilience
 * @param w     width in 8 pixel blocks
 * @param h     height in 8 pixel blocks
 */
static void v_block_filter(ERContext *s, uint8_t *dst, int w, int h,
                           ptrdiff_t stride, int is_luma)
{
    int b_x, b_y;
    ptrdiff_t mvx_stride, mvy_stride;
    const uint8_t *cm = ff_crop_tab + MAX_NEG_CROP;
    set_mv_strides(s, &mvx_stride, &mvy_stride);
    mvx_stride >>= is_luma;
    mvy_stride *= mvx_stride;

    for (b_y = 0; b_y < h - 1; b_y++) {
        for (b_x = 0; b_x < w; b_x++) {
            int x;
            int top_status    = s->error_status_table[(b_x >> is_luma) +  (b_y      >> is_luma) * s->mb_stride];
            int bottom_status = s->error_status_table[(b_x >> is_luma) + ((b_y + 1) >> is_luma) * s->mb_stride];
            int top_intra     = IS_INTRA(s->cur_pic.mb_type[(b_x >> is_luma) + ( b_y      >> is_luma) * s->mb_stride]);
            int bottom_intra  = IS_INTRA(s->cur_pic.mb_type[(b_x >> is_luma) + ((b_y + 1) >> is_luma) * s->mb_stride]);
            int top_damage    = top_status & ER_MB_ERROR;
            int bottom_damage = bottom_status & ER_MB_ERROR;
            int offset        = b_x * 8 + b_y * stride * 8;

            int16_t *top_mv    = s->cur_pic.motion_val[0][mvy_stride *  b_y      + mvx_stride * b_x];
            int16_t *bottom_mv = s->cur_pic.motion_val[0][mvy_stride * (b_y + 1) + mvx_stride * b_x];

            if (!(top_damage || bottom_damage))
                continue; // both undamaged

            if ((!top_intra) && (!bottom_intra) &&
                FFABS(top_mv[0] - bottom_mv[0]) +
                FFABS(top_mv[1] + bottom_mv[1]) < 2)
                continue;

            for (x = 0; x < 8; x++) {
                int a, b, c, d;

                a = dst[offset + x + 7 * stride] - dst[offset + x + 6 * stride];
                b = dst[offset + x + 8 * stride] - dst[offset + x + 7 * stride];
                c = dst[offset + x + 9 * stride] - dst[offset + x + 8 * stride];

                d = FFABS(b) - ((FFABS(a) + FFABS(c) + 1) >> 1);
                d = FFMAX(d, 0);
                if (b < 0)
                    d = -d;

                if (d == 0)
                    continue;

                if (!(top_damage && bottom_damage))
                    d = d * 16 / 9;

                if (top_damage) {
                    dst[offset + x +  7 * stride] = cm[dst[offset + x +  7 * stride] + ((d * 7) >> 4)];
                    dst[offset + x +  6 * stride] = cm[dst[offset + x +  6 * stride] + ((d * 5) >> 4)];
                    dst[offset + x +  5 * stride] = cm[dst[offset + x +  5 * stride] + ((d * 3) >> 4)];
                    dst[offset + x +  4 * stride] = cm[dst[offset + x +  4 * stride] + ((d * 1) >> 4)];
                }
                if (bottom_damage) {
                    dst[offset + x +  8 * stride] = cm[dst[offset + x +  8 * stride] - ((d * 7) >> 4)];
                    dst[offset + x +  9 * stride] = cm[dst[offset + x +  9 * stride] - ((d * 5) >> 4)];
                    dst[offset + x + 10 * stride] = cm[dst[offset + x + 10 * stride] - ((d * 3) >> 4)];
                    dst[offset + x + 11 * stride] = cm[dst[offset + x + 11 * stride] - ((d * 1) >> 4)];
                }
            }
        }
    }
}

#define MV_FROZEN    8
#define MV_CHANGED   4
#define MV_UNCHANGED 2
#define MV_LISTED    1
static av_always_inline void add_blocklist(int (*blocklist)[2], int *blocklist_length, uint8_t *fixed, int mb_x, int mb_y, int mb_xy)
{
    if (fixed[mb_xy])
        return;
    fixed[mb_xy] = MV_LISTED;
    blocklist[ *blocklist_length   ][0] = mb_x;
    blocklist[(*blocklist_length)++][1] = mb_y;
}

static void guess_mv(ERContext *s)
{
    int (*blocklist)[2], (*next_blocklist)[2];
    uint8_t *fixed;
    const ptrdiff_t mb_stride = s->mb_stride;
    const int mb_width  = s->mb_width;
    int mb_height = s->mb_height;
    int i, depth, num_avail;
    int mb_x, mb_y;
    ptrdiff_t mot_step, mot_stride;
    int blocklist_length, next_blocklist_length;

    if (s->last_pic.f && s->last_pic.f->data[0])
        mb_height = FFMIN(mb_height, (s->last_pic.f->height+15)>>4);
    if (s->next_pic.f && s->next_pic.f->data[0])
        mb_height = FFMIN(mb_height, (s->next_pic.f->height+15)>>4);

    blocklist      = (int (*)[2])s->er_temp_buffer;
    next_blocklist = blocklist + s->mb_stride * s->mb_height;
    fixed          = (uint8_t *)(next_blocklist + s->mb_stride * s->mb_height);

    set_mv_strides(s, &mot_step, &mot_stride);

    num_avail = 0;
    if (s->last_pic.motion_val[0])
        ff_thread_await_progress(s->last_pic.tf, mb_height-1, 0);
    for (i = 0; i < mb_width * mb_height; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int f = 0;
        int error = s->error_status_table[mb_xy];

        if (IS_INTRA(s->cur_pic.mb_type[mb_xy]))
            f = MV_FROZEN; // intra // FIXME check
        if (!(error & ER_MV_ERROR))
            f = MV_FROZEN; // inter with undamaged MV

        fixed[mb_xy] = f;
        if (f == MV_FROZEN)
            num_avail++;
        else if(s->last_pic.f->data[0] && s->last_pic.motion_val[0]){
            const int mb_y= mb_xy / s->mb_stride;
            const int mb_x= mb_xy % s->mb_stride;
            const int mot_index= (mb_x + mb_y*mot_stride) * mot_step;
            s->cur_pic.motion_val[0][mot_index][0]= s->last_pic.motion_val[0][mot_index][0];
            s->cur_pic.motion_val[0][mot_index][1]= s->last_pic.motion_val[0][mot_index][1];
            s->cur_pic.ref_index[0][4*mb_xy]      = s->last_pic.ref_index[0][4*mb_xy];
        }
    }

    if ((!(s->avctx->error_concealment&FF_EC_GUESS_MVS)) ||
        num_avail <= FFMAX(mb_width, mb_height) / 2) {
        for (mb_y = 0; mb_y < mb_height; mb_y++) {
            for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
                const int mb_xy = mb_x + mb_y * s->mb_stride;
                int mv_dir = (s->last_pic.f && s->last_pic.f->data[0]) ? MV_DIR_FORWARD : MV_DIR_BACKWARD;

                if (IS_INTRA(s->cur_pic.mb_type[mb_xy]))
                    continue;
                if (!(s->error_status_table[mb_xy] & ER_MV_ERROR))
                    continue;

                s->mv[0][0][0] = 0;
                s->mv[0][0][1] = 0;
                s->decode_mb(s->opaque, 0, mv_dir, MV_TYPE_16X16, &s->mv,
                             mb_x, mb_y, 0, 0);
            }
        }
        return;
    }

    blocklist_length = 0;
    for (mb_y = 0; mb_y < mb_height; mb_y++) {
        for (mb_x = 0; mb_x < mb_width; mb_x++) {
            const int mb_xy = mb_x + mb_y * mb_stride;
            if (fixed[mb_xy] == MV_FROZEN) {
                if (mb_x)               add_blocklist(blocklist, &blocklist_length, fixed, mb_x - 1, mb_y, mb_xy - 1);
                if (mb_y)               add_blocklist(blocklist, &blocklist_length, fixed, mb_x, mb_y - 1, mb_xy - mb_stride);
                if (mb_x+1 < mb_width)  add_blocklist(blocklist, &blocklist_length, fixed, mb_x + 1, mb_y, mb_xy + 1);
                if (mb_y+1 < mb_height) add_blocklist(blocklist, &blocklist_length, fixed, mb_x, mb_y + 1, mb_xy + mb_stride);
            }
        }
    }

    for (depth = 0; ; depth++) {
        int changed, pass, none_left;
        int blocklist_index;

        none_left = 1;
        changed   = 1;
        for (pass = 0; (changed || pass < 2) && pass < 10; pass++) {
            changed = 0;
            for (blocklist_index = 0; blocklist_index < blocklist_length; blocklist_index++) {
                const int mb_x = blocklist[blocklist_index][0];
                const int mb_y = blocklist[blocklist_index][1];
                const int mb_xy = mb_x + mb_y * mb_stride;
                int mv_predictor[8][2];
                int ref[8];
                int pred_count;
                int j;
                int best_score;
                int best_pred;
                int mot_index;
                int prev_x, prev_y, prev_ref;

                if ((mb_x ^ mb_y ^ pass) & 1)
                    continue;
                av_assert2(fixed[mb_xy] != MV_FROZEN);


                av_assert1(!IS_INTRA(s->cur_pic.mb_type[mb_xy]));
                av_assert1(s->last_pic.f && s->last_pic.f->data[0]);

                j = 0;
                if (mb_x > 0)
                    j |= fixed[mb_xy - 1];
                if (mb_x + 1 < mb_width)
                    j |= fixed[mb_xy + 1];
                if (mb_y > 0)
                    j |= fixed[mb_xy - mb_stride];
                if (mb_y + 1 < mb_height)
                    j |= fixed[mb_xy + mb_stride];

                av_assert2(j & MV_FROZEN);

                if (!(j & MV_CHANGED) && pass > 1)
                    continue;

                none_left = 0;
                pred_count = 0;
                mot_index  = (mb_x + mb_y * mot_stride) * mot_step;

                if (mb_x > 0 && fixed[mb_xy - 1] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index - mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index - mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy - 1)];
                    pred_count++;
                }
                if (mb_x + 1 < mb_width && fixed[mb_xy + 1] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index + mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index + mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy + 1)];
                    pred_count++;
                }
                if (mb_y > 0 && fixed[mb_xy - mb_stride] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index - mot_stride * mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index - mot_stride * mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy - s->mb_stride)];
                    pred_count++;
                }
                if (mb_y + 1<mb_height && fixed[mb_xy + mb_stride] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index + mot_stride * mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index + mot_stride * mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy + s->mb_stride)];
                    pred_count++;
                }
                if (pred_count == 0)
                    continue;

                if (pred_count > 1) {
                    int sum_x = 0, sum_y = 0, sum_r = 0;
                    int max_x, max_y, min_x, min_y, max_r, min_r;

                    for (j = 0; j < pred_count; j++) {
                        sum_x += mv_predictor[j][0];
                        sum_y += mv_predictor[j][1];
                        sum_r += ref[j];
                        if (j && ref[j] != ref[j - 1])
                            goto skip_mean_and_median;
                    }

                    /* mean */
                    mv_predictor[pred_count][0] = sum_x / j;
                    mv_predictor[pred_count][1] = sum_y / j;
                             ref[pred_count]    = sum_r / j;

                    /* median */
                    if (pred_count >= 3) {
                        min_y = min_x = min_r =  99999;
                        max_y = max_x = max_r = -99999;
                    } else {
                        min_x = min_y = max_x = max_y = min_r = max_r = 0;
                    }
                    for (j = 0; j < pred_count; j++) { 
                        max_x = FFMAX(max_x, mv_predictor[j][0]);
                        max_y = FFMAX(max_y, mv_predictor[j][1]);
                        max_r = FFMAX(max_r, ref[j]);
                        min_x = FFMIN(min_x, mv_predictor[j][0]);
                        min_y = FFMIN(min_y, mv_predictor[j][1]);
                        min_r = FFMIN(min_r, ref[j]);
                    }
                    mv_predictor[pred_count + 1][0] = sum_x - max_x - min_x;
                    mv_predictor[pred_count + 1][1] = sum_y - max_y - min_y;
                             ref[pred_count + 1]    = sum_r - max_r - min_r;

                    if (pred_count == 4) {
                        mv_predictor[pred_count + 1][0] /= 2;
                        mv_predictor[pred_count + 1][1] /= 2;
                                 ref[pred_count + 1]    /= 2;
                    }
                    pred_count += 2;
                }

skip_mean_and_median:
                /* zero MV */
                mv_predictor[pred_count][0] =
                mv_predictor[pred_count][1] =
                         ref[pred_count]    = 0;
                pred_count++;

                prev_x   = s->cur_pic.motion_val[0][mot_index][0];
                prev_y   = s->cur_pic.motion_val[0][mot_index][1];
                prev_ref = s->cur_pic.ref_index[0][4 * mb_xy];

                /* last MV */
                mv_predictor[pred_count][0] = prev_x;
                mv_predictor[pred_count][1] = prev_y;
                         ref[pred_count]    = prev_ref;
                pred_count++;

                best_pred = 0;
                best_score = 256 * 256 * 256 * 64;
                for (j = 0; j < pred_count; j++) {
                    int *linesize = s->cur_pic.f->linesize;
                    int score = 0;
                    uint8_t *src = s->cur_pic.f->data[0] +
                                   mb_x * 16 + mb_y * 16 * linesize[0];

                    s->cur_pic.motion_val[0][mot_index][0] =
                        s->mv[0][0][0] = mv_predictor[j][0];
                    s->cur_pic.motion_val[0][mot_index][1] =
                        s->mv[0][0][1] = mv_predictor[j][1];

                    // predictor intra or otherwise not available
                    if (ref[j] < 0)
                        continue;

                    s->decode_mb(s->opaque, ref[j], MV_DIR_FORWARD,
                                 MV_TYPE_16X16, &s->mv, mb_x, mb_y, 0, 0);

                    if (mb_x > 0 && fixed[mb_xy - 1] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(src[k * linesize[0] - 1] -
                                           src[k * linesize[0]]);
                    }
                    if (mb_x + 1 < mb_width && fixed[mb_xy + 1] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(src[k * linesize[0] + 15] -
                                           src[k * linesize[0] + 16]);
                    }
                    if (mb_y > 0 && fixed[mb_xy - mb_stride] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(src[k - linesize[0]] - src[k]);
                    }
                    if (mb_y + 1 < mb_height && fixed[mb_xy + mb_stride] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(src[k + linesize[0] * 15] -
                                           src[k + linesize[0] * 16]);
                    }

                    if (score <= best_score) { // <= will favor the last MV
                        best_score = score;
                        best_pred  = j;
                    }
                }
                s->mv[0][0][0] = mv_predictor[best_pred][0];
                s->mv[0][0][1] = mv_predictor[best_pred][1];

                for (i = 0; i < mot_step; i++)
                    for (j = 0; j < mot_step; j++) {
                        s->cur_pic.motion_val[0][mot_index + i + j * mot_stride][0] = s->mv[0][0][0];
                        s->cur_pic.motion_val[0][mot_index + i + j * mot_stride][1] = s->mv[0][0][1];
                    }

                s->decode_mb(s->opaque, ref[best_pred], MV_DIR_FORWARD,
                             MV_TYPE_16X16, &s->mv, mb_x, mb_y, 0, 0);


                if (s->mv[0][0][0] != prev_x || s->mv[0][0][1] != prev_y) {
                    fixed[mb_xy] = MV_CHANGED;
                    changed++;
                } else
                    fixed[mb_xy] = MV_UNCHANGED;
            }
        }

        if (none_left)
            return;

        next_blocklist_length = 0;

        for (blocklist_index = 0; blocklist_index < blocklist_length; blocklist_index++) {
            const int mb_x = blocklist[blocklist_index][0];
            const int mb_y = blocklist[blocklist_index][1];
            const int mb_xy = mb_x + mb_y * mb_stride;

            if (fixed[mb_xy] & (MV_CHANGED|MV_UNCHANGED|MV_FROZEN)) {
                fixed[mb_xy] = MV_FROZEN;
                if (mb_x > 0)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x - 1, mb_y, mb_xy - 1);
                if (mb_y > 0)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x, mb_y - 1, mb_xy - mb_stride);
                if (mb_x + 1 < mb_width)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x + 1, mb_y, mb_xy + 1);
                if (mb_y + 1 < mb_height)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x, mb_y + 1, mb_xy + mb_stride);
            }
        }
        av_assert0(next_blocklist_length <= mb_height * mb_width);
        FFSWAP(int , blocklist_length, next_blocklist_length);
        FFSWAP(void*, blocklist, next_blocklist);
    }
}

static void guess_mv_edge(ERContext *s)
{


    
    //s->cur_pic.edge_mag = av_malloc_array(s->cur_pic.f->width, s->cur_pic.f->height*sizeof(float));
   // s->cur_pic.edge_dir = av_malloc_array(s->cur_pic.f->width, s->cur_pic.f->height*sizeof(float));
   
    int (*blocklist)[2], (*next_blocklist)[2];
    uint8_t *fixed;
    const ptrdiff_t mb_stride = s->mb_stride;
    const int mb_width  = s->mb_width;
    int mb_height = s->mb_height;
    int i, depth, num_avail;
    int mb_x, mb_y;
    ptrdiff_t mot_step, mot_stride;
    int blocklist_length, next_blocklist_length;

    if (s->last_pic.f && s->last_pic.f->data[0])
        mb_height = FFMIN(mb_height, (s->last_pic.f->height+15)>>4);
    if (s->next_pic.f && s->next_pic.f->data[0])
        mb_height = FFMIN(mb_height, (s->next_pic.f->height+15)>>4);

    blocklist      = (int (*)[2])s->er_temp_buffer;
    next_blocklist = blocklist + s->mb_stride * s->mb_height;
    fixed          = (uint8_t *)(next_blocklist + s->mb_stride * s->mb_height);

    set_mv_strides(s, &mot_step, &mot_stride);

    num_avail = 0;
    if (s->last_pic.motion_val[0])
        ff_thread_await_progress(s->last_pic.tf, mb_height-1, 0);

    // kjlee: mark normal macroblocks (intra and no mv error) as frozen
    // kjlee: for mb with mv error, if last pic mv is normal copy
    for (i = 0; i < mb_width * mb_height; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int f = 0;
        int error = s->error_status_table[mb_xy];

        if (IS_INTRA(s->cur_pic.mb_type[mb_xy]))
            f = MV_FROZEN; // intra // FIXME check
        if (!(error & ER_MV_ERROR))
            f = MV_FROZEN; // inter with undamaged MV

        fixed[mb_xy] = f;
        if (f == MV_FROZEN)
            num_avail++;
        else if(s->last_pic.f->data[0] && s->last_pic.motion_val[0]){
            const int mb_y= mb_xy / s->mb_stride;
            const int mb_x= mb_xy % s->mb_stride;
            const int mot_index= (mb_x + mb_y*mot_stride) * mot_step;
            s->cur_pic.motion_val[0][mot_index][0]= s->last_pic.motion_val[0][mot_index][0];
            s->cur_pic.motion_val[0][mot_index][1]= s->last_pic.motion_val[0][mot_index][1];
            s->cur_pic.ref_index[0][4*mb_xy]      = s->last_pic.ref_index[0][4*mb_xy];
        }
    }


    // kjlee: if error conceal with guess mv turned off, set mv as zeros and return
    if ((!(s->avctx->error_concealment&FF_EC_GUESS_MVS)) ||
        num_avail <= FFMAX(mb_width, mb_height) / 2) {
        for (mb_y = 0; mb_y < mb_height; mb_y++) {
            for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
                const int mb_xy = mb_x + mb_y * s->mb_stride;
                int mv_dir = (s->last_pic.f && s->last_pic.f->data[0]) ? MV_DIR_FORWARD : MV_DIR_BACKWARD;

                if (IS_INTRA(s->cur_pic.mb_type[mb_xy]))
                    continue;
                if (!(s->error_status_table[mb_xy] & ER_MV_ERROR))
                    continue;

                s->mv[0][0][0] = 0;
                s->mv[0][0][1] = 0;
                s->decode_mb(s->opaque, 0, mv_dir, MV_TYPE_16X16, &s->mv,
                             mb_x, mb_y, 0, 0);
            }
        }
        return;
    }

    blocklist_length = 0;
    for (mb_y = 0; mb_y < mb_height; mb_y++) {
        for (mb_x = 0; mb_x < mb_width; mb_x++) {
            const int mb_xy = mb_x + mb_y * mb_stride;
            if (fixed[mb_xy] == MV_FROZEN) {
                if (mb_x)               add_blocklist(blocklist, &blocklist_length, fixed, mb_x - 1, mb_y, mb_xy - 1);
                if (mb_y)               add_blocklist(blocklist, &blocklist_length, fixed, mb_x, mb_y - 1, mb_xy - mb_stride);
                if (mb_x+1 < mb_width)  add_blocklist(blocklist, &blocklist_length, fixed, mb_x + 1, mb_y, mb_xy + 1);
                if (mb_y+1 < mb_height) add_blocklist(blocklist, &blocklist_length, fixed, mb_x, mb_y + 1, mb_xy + mb_stride);
            }
        }
    }

    for (depth = 0; ; depth++) {
        int changed, pass, none_left;
        int blocklist_index;

        none_left = 1;
        changed   = 1;
        for (pass = 0; (changed || pass < 2) && pass < 10; pass++) {
            changed = 0;
            for (blocklist_index = 0; blocklist_index < blocklist_length; blocklist_index++) {
                const int mb_x = blocklist[blocklist_index][0];
                const int mb_y = blocklist[blocklist_index][1];
                const int mb_xy = mb_x + mb_y * mb_stride;
                int mv_predictor[8][2];
                int ref[8];
                int pred_count;
                int j;
                int best_score;
                int best_pred;
                int mot_index;
                int prev_x, prev_y, prev_ref;

                if ((mb_x ^ mb_y ^ pass) & 1)
                    continue;
                av_assert2(fixed[mb_xy] != MV_FROZEN);


                av_assert1(!IS_INTRA(s->cur_pic.mb_type[mb_xy]));
                av_assert1(s->last_pic.f && s->last_pic.f->data[0]);

                j = 0;
                if (mb_x > 0)
                    j |= fixed[mb_xy - 1];
                if (mb_x + 1 < mb_width)
                    j |= fixed[mb_xy + 1];
                if (mb_y > 0)
                    j |= fixed[mb_xy - mb_stride];
                if (mb_y + 1 < mb_height)
                    j |= fixed[mb_xy + mb_stride];

                av_assert2(j & MV_FROZEN);

                if (!(j & MV_CHANGED) && pass > 1)
                    continue;

                none_left = 0;
                pred_count = 0;
                mot_index  = (mb_x + mb_y * mot_stride) * mot_step;

                // left mb
                if (mb_x > 0 && fixed[mb_xy - 1] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index - mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index - mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy - 1)];
                    pred_count++;
                }
                // right mb
                if (mb_x + 1 < mb_width && fixed[mb_xy + 1] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index + mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index + mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy + 1)];
                    pred_count++;
                }
                // top mb
                if (mb_y > 0 && fixed[mb_xy - mb_stride] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index - mot_stride * mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index - mot_stride * mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy - s->mb_stride)];
                    pred_count++;
                }
                // bottom mb
                if (mb_y + 1<mb_height && fixed[mb_xy + mb_stride] > 1) {
                    mv_predictor[pred_count][0] =
                        s->cur_pic.motion_val[0][mot_index + mot_stride * mot_step][0];
                    mv_predictor[pred_count][1] =
                        s->cur_pic.motion_val[0][mot_index + mot_stride * mot_step][1];
                    ref[pred_count] =
                        s->cur_pic.ref_index[0][4 * (mb_xy + s->mb_stride)];
                    pred_count++;
                }
                if (pred_count == 0)
                    continue;

                if (pred_count > 1) {
                    int sum_x = 0, sum_y = 0, sum_r = 0;
                    int max_x, max_y, min_x, min_y, max_r, min_r;

                    for (j = 0; j < pred_count; j++) {
                        sum_x += mv_predictor[j][0];
                        sum_y += mv_predictor[j][1];
                        sum_r += ref[j];
                        if (j && ref[j] != ref[j - 1])
                            goto skip_mean_and_median;
                    }

                    /* mean */
                    mv_predictor[pred_count][0] = sum_x / j;
                    mv_predictor[pred_count][1] = sum_y / j;
                             ref[pred_count]    = sum_r / j;

                    /* median */
                    if (pred_count >= 3) {
                        min_y = min_x = min_r =  99999;
                        max_y = max_x = max_r = -99999;
                    } else {
                        min_x = min_y = max_x = max_y = min_r = max_r = 0;
                    }
                    for (j = 0; j < pred_count; j++) {
                        max_x = FFMAX(max_x, mv_predictor[j][0]);
                        max_y = FFMAX(max_y, mv_predictor[j][1]);
                        max_r = FFMAX(max_r, ref[j]);
                        min_x = FFMIN(min_x, mv_predictor[j][0]);
                        min_y = FFMIN(min_y, mv_predictor[j][1]);
                        min_r = FFMIN(min_r, ref[j]);
                    }
                    mv_predictor[pred_count + 1][0] = sum_x - max_x - min_x;
                    mv_predictor[pred_count + 1][1] = sum_y - max_y - min_y;
                             ref[pred_count + 1]    = sum_r - max_r - min_r;

                    if (pred_count == 4) {
                        mv_predictor[pred_count + 1][0] /= 2;
                        mv_predictor[pred_count + 1][1] /= 2;
                                 ref[pred_count + 1]    /= 2;
                    }
                    pred_count += 2;
                }

skip_mean_and_median:
                /* zero MV */
                mv_predictor[pred_count][0] =
                mv_predictor[pred_count][1] =
                         ref[pred_count]    = 0;
                pred_count++;

                prev_x   = s->cur_pic.motion_val[0][mot_index][0];
                prev_y   = s->cur_pic.motion_val[0][mot_index][1];
                prev_ref = s->cur_pic.ref_index[0][4 * mb_xy];

                /* last MV */
                mv_predictor[pred_count][0] = prev_x;
                mv_predictor[pred_count][1] = prev_y;
                         ref[pred_count]    = prev_ref;
                pred_count++;

                // kjlee: TODO: check edge dominance of nearby mbs

                best_pred = 0;
                best_score = 256 * 256 * 256 * 64;
                for (j = 0; j < pred_count; j++) {
                    int *linesize = s->cur_pic.f->linesize;
                    int score = 0;
                    uint8_t *src = s->cur_pic.f->data[0] +
                                   mb_x * 16 + mb_y * 16 * linesize[0];

                    s->cur_pic.motion_val[0][mot_index][0] =
                        s->mv[0][0][0] = mv_predictor[j][0];
                    s->cur_pic.motion_val[0][mot_index][1] =
                        s->mv[0][0][1] = mv_predictor[j][1];

                    // predictor intra or otherwise not available
                    if (ref[j] < 0)
                        continue;

                    s->decode_mb(s->opaque, ref[j], MV_DIR_FORWARD,
                                 MV_TYPE_16X16, &s->mv, mb_x, mb_y, 0, 0);

                    // kjlee: update edge detection for current macroblock with yuv value of new estimated motion vector
                    edge_detection_mb(s, mb_x, mb_y); 
                    float *edge = s->cur_pic.edge_dir +
                                   mb_x * 16 + mb_y * 16 * linesize[0];
                    // kjlee: currently simply calculate edge direction difference
                    if (mb_x > 0 && fixed[mb_xy - 1] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(edge[k * linesize[0] - 1] -
                                           edge[k * linesize[0]]);
                    }
                    if (mb_x + 1 < mb_width && fixed[mb_xy + 1] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(edge[k * linesize[0] + 15] -
                                           edge[k * linesize[0] + 16]);
                    }
                    if (mb_y > 0 && fixed[mb_xy - mb_stride] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(edge[k - linesize[0]] - edge[k]);
                    }
                    if (mb_y + 1 < mb_height && fixed[mb_xy + mb_stride] > 1) {
                        int k;
                        for (k = 0; k < 16; k++)
                            score += FFABS(edge[k + linesize[0] * 15] -
                                           edge[k + linesize[0] * 16]);
                    }

                    if (score <= best_score) { // <= will favor the last MV
                        best_score = score;
                        best_pred  = j;
                    }
                }
                s->mv[0][0][0] = mv_predictor[best_pred][0];
                s->mv[0][0][1] = mv_predictor[best_pred][1];

                for (i = 0; i < mot_step; i++)
                    for (j = 0; j < mot_step; j++) {
                        s->cur_pic.motion_val[0][mot_index + i + j * mot_stride][0] = s->mv[0][0][0];
                        s->cur_pic.motion_val[0][mot_index + i + j * mot_stride][1] = s->mv[0][0][1];
                    }

                s->decode_mb(s->opaque, ref[best_pred], MV_DIR_FORWARD,
                             MV_TYPE_16X16, &s->mv, mb_x, mb_y, 0, 0);


                if (s->mv[0][0][0] != prev_x || s->mv[0][0][1] != prev_y) {
                    fixed[mb_xy] = MV_CHANGED;
                    changed++;
                } else
                    fixed[mb_xy] = MV_UNCHANGED;
            }
        }

        if (none_left)
            return;

        next_blocklist_length = 0;

        for (blocklist_index = 0; blocklist_index < blocklist_length; blocklist_index++) {
            const int mb_x = blocklist[blocklist_index][0];
            const int mb_y = blocklist[blocklist_index][1];
            const int mb_xy = mb_x + mb_y * mb_stride;

            if (fixed[mb_xy] & (MV_CHANGED|MV_UNCHANGED|MV_FROZEN)) {
                fixed[mb_xy] = MV_FROZEN;
                if (mb_x > 0)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x - 1, mb_y, mb_xy - 1);
                if (mb_y > 0)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x, mb_y - 1, mb_xy - mb_stride);
                if (mb_x + 1 < mb_width)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x + 1, mb_y, mb_xy + 1);
                if (mb_y + 1 < mb_height)
                    add_blocklist(next_blocklist, &next_blocklist_length, fixed, mb_x, mb_y + 1, mb_xy + mb_stride);
            }
        }
        av_assert0(next_blocklist_length <= mb_height * mb_width);
        FFSWAP(int , blocklist_length, next_blocklist_length);
        FFSWAP(void*, blocklist, next_blocklist);
    }
}


static int is_intra_more_likely(ERContext *s)
{
    int is_intra_likely, i, j, undamaged_count, skip_amount, mb_x, mb_y;

    if (!s->last_pic.f || !s->last_pic.f->data[0])
        return 1; // no previous frame available -> use spatial prediction

    if (s->avctx->error_concealment & FF_EC_FAVOR_INTER)
        return 0;

    undamaged_count = 0;
    for (i = 0; i < s->mb_num; i++) {
        const int mb_xy = s->mb_index2xy[i];
        const int error = s->error_status_table[mb_xy];
        if (!((error & ER_DC_ERROR) && (error & ER_MV_ERROR)))
            undamaged_count++;
    }

    if (undamaged_count < 5)
        return 0; // almost all MBs damaged -> use temporal prediction

    skip_amount     = FFMAX(undamaged_count / 50, 1); // check only up to 50 MBs
    is_intra_likely = 0;

    j = 0;
    for (mb_y = 0; mb_y < s->mb_height - 1; mb_y++) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
            int error;
            const int mb_xy = mb_x + mb_y * s->mb_stride;

            error = s->error_status_table[mb_xy];
            if ((error & ER_DC_ERROR) && (error & ER_MV_ERROR))
                continue; // skip damaged

            j++;
            // skip a few to speed things up
            if ((j % skip_amount) != 0)
                continue;

            if (s->cur_pic.f->pict_type == AV_PICTURE_TYPE_I) {
                int *linesize = s->cur_pic.f->linesize;
                uint8_t *mb_ptr      = s->cur_pic.f->data[0] +
                                       mb_x * 16 + mb_y * 16 * linesize[0];
                uint8_t *last_mb_ptr = s->last_pic.f->data[0] +
                                       mb_x * 16 + mb_y * 16 * linesize[0];

                if (s->avctx->codec_id == AV_CODEC_ID_H264) {
                    // FIXME
                } else {
                    ff_thread_await_progress(s->last_pic.tf, mb_y, 0);
                }
                is_intra_likely += s->mecc.sad[0](NULL, last_mb_ptr, mb_ptr,
                                                  linesize[0], 16);
                // FIXME need await_progress() here
                is_intra_likely -= s->mecc.sad[0](NULL, last_mb_ptr,
                                                  last_mb_ptr + linesize[0] * 16,
                                                  linesize[0], 16);
            } else {
                if (IS_INTRA(s->cur_pic.mb_type[mb_xy]))
                   is_intra_likely++;
                else
                   is_intra_likely--;
            }
        }
    }
//      av_log(NULL, AV_LOG_ERROR, "is_intra_likely: %d type:%d\n", is_intra_likely, s->pict_type);
    return is_intra_likely > 0;
}


// initialize error_status_table
void ff_er_frame_start(ERContext *s)
{
    if (!s->avctx->error_concealment)
        return;

    if (!s->mecc_inited) {
        ff_me_cmp_init(&s->mecc, s->avctx);
        s->mecc_inited = 1;
    }

    // ER_MB_ERROR == (ER_AC_ERROR|ER_DC_ERROR|ER_MV_ERROR) == (2 | 4 | 8)  == 0000 1110
    // VP_START  == 1 == 0000 0001
    // ER_MB_END == (ER_AC_END|ER_DC_END|ER_MV_END) == (16 | 32 | 64) == 0111 000
    // ER_MB_ERROR | VP_START | ER_MB_END = 0111 1111
    memset(s->error_status_table, ER_MB_ERROR | VP_START | ER_MB_END,
           s->mb_stride * s->mb_height * sizeof(uint8_t));
    atomic_init(&s->error_count, 3 * s->mb_num);
    s->error_occurred = 0;
}

static int er_supported(ERContext *s)
{
    if(s->avctx->hwaccel && s->avctx->hwaccel->decode_slice           ||
       !s->cur_pic.f                                                  ||
       s->cur_pic.field_picture
    )
        return 0;
    return 1;
}

/**
 * Add a slice.
 * @param endx   x component of the last macroblock, can be -1
 *               for the last of the previous line
 * @param status the status at the end (ER_MV_END, ER_AC_ERROR, ...), it is
 *               assumed that no earlier end or error of the same type occurred
 */


int temp = 0;
void ff_er_add_slice(ERContext *s, int startx, int starty,
                     int endx, int endy, int status)
{
    const int start_i  = av_clip(startx + starty * s->mb_width, 0, s->mb_num - 1);
    const int end_i    = av_clip(endx   + endy   * s->mb_width, 0, s->mb_num);
    const int start_xy = s->mb_index2xy[start_i];
    const int end_xy   = s->mb_index2xy[end_i];
    int mask           = -1;

    if (s->avctx->hwaccel && s->avctx->hwaccel->decode_slice)
        return;

    if (start_i > end_i || start_xy > end_xy) {
        av_log(s->avctx, AV_LOG_ERROR,
               "internal error, slice end before start\n");
        return;
    }

    if (!s->avctx->error_concealment)
        return;

    mask &= ~VP_START;
    if (status & (ER_AC_ERROR | ER_AC_END)) {
        mask           &= ~(ER_AC_ERROR | ER_AC_END);
        atomic_fetch_add(&s->error_count, start_i - end_i - 1);
    }
    if (status & (ER_DC_ERROR | ER_DC_END)) {
        mask           &= ~(ER_DC_ERROR | ER_DC_END);
        atomic_fetch_add(&s->error_count, start_i - end_i - 1);
    }
    if (status & (ER_MV_ERROR | ER_MV_END)) {
        mask           &= ~(ER_MV_ERROR | ER_MV_END);
        atomic_fetch_add(&s->error_count, start_i - end_i - 1);
    }

    if (status & ER_MB_ERROR) {
        s->error_occurred = 1;
        atomic_store(&s->error_count, INT_MAX);
    }

    if (mask == ~0x7F) {
        memset(&s->error_status_table[start_xy], 0,
               (end_xy - start_xy) * sizeof(uint8_t));
    } else {
        int i;
        for (i = start_xy; i < end_xy; i++)
            s->error_status_table[i] &= mask;
    }

    if (end_i == s->mb_num)
        atomic_store(&s->error_count, INT_MAX);
    else {
        s->error_status_table[end_xy] &= mask;
        s->error_status_table[end_xy] |= status;
    }

    s->error_status_table[start_xy] |= VP_START;

    if (start_xy > 0 && !(s->avctx->active_thread_type & FF_THREAD_SLICE) &&
        er_supported(s) && s->avctx->skip_top * s->mb_width < start_i) {
        int prev_status = s->error_status_table[s->mb_index2xy[start_i - 1]];

        prev_status &= ~ VP_START;
        if (prev_status != (ER_MV_END | ER_DC_END | ER_AC_END)) {
            s->error_occurred = 1;
            atomic_store(&s->error_count, INT_MAX);
        }
    }
}


// this function called, performs error concealment
void ff_er_frame_end(ERContext *s)
{
    int *linesize = NULL;
    int i, mb_x, mb_y, error, error_type, dc_error, mv_error, ac_error;
    int distance;
    int threshold_part[4] = { 100, 100, 100 };
    int threshold = 50;
    int is_intra_likely;
    int size = s->b8_stride * 2 * s->mb_height;

    /* We do not support ER of field pictures yet,
     * though it should not crash if enabled. */
    if (!s->avctx->error_concealment || !atomic_load(&s->error_count)  ||
        s->avctx->lowres                                               ||
        !er_supported(s)                                               ||
        atomic_load(&s->error_count) == 3 * s->mb_width *
                          (s->avctx->skip_top + s->avctx->skip_bottom)) {
        return;
    }
    linesize = s->cur_pic.f->linesize;

    if (   s->avctx->codec_id == AV_CODEC_ID_MPEG2VIDEO
        && (FFALIGN(s->avctx->height, 16)&16)
        && atomic_load(&s->error_count) == 3 * s->mb_width * (s->avctx->skip_top + s->avctx->skip_bottom + 1)) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
            int status = s->error_status_table[mb_x + (s->mb_height - 1) * s->mb_stride];
            if (status != 0x7F)
                break;
        }

        if (mb_x == s->mb_width) {
            av_log(s->avctx, AV_LOG_DEBUG, "ignoring last missing slice\n");
            return;
        }
    }

    // if format/size of last picture is not the same as current picture
    if (s->last_pic.f) {
        if (s->last_pic.f->width  != s->cur_pic.f->width  ||
            s->last_pic.f->height != s->cur_pic.f->height ||
            s->last_pic.f->format != s->cur_pic.f->format) {
            av_log(s->avctx, AV_LOG_WARNING, "Cannot use previous picture in error concealment\n");
            memset(&s->last_pic, 0, sizeof(s->last_pic));
        }
    }
    // if format/size of next picture is not the same as current picture
    if (s->next_pic.f) {
        if (s->next_pic.f->width  != s->cur_pic.f->width  ||
            s->next_pic.f->height != s->cur_pic.f->height ||
            s->next_pic.f->format != s->cur_pic.f->format) {
            av_log(s->avctx, AV_LOG_WARNING, "Cannot use next picture in error concealment\n");
            memset(&s->next_pic, 0, sizeof(s->next_pic));
        }
    }

    if (!s->cur_pic.motion_val[0] || !s->cur_pic.ref_index[0]) {
        av_log(s->avctx, AV_LOG_ERROR, "Warning MVs not available\n");

        for (i = 0; i < 2; i++) {
            s->ref_index_buf[i]  = av_buffer_allocz(s->mb_stride * s->mb_height * 4 * sizeof(uint8_t));
            s->motion_val_buf[i] = av_buffer_allocz((size + 4) * 2 * sizeof(uint16_t));
            if (!s->ref_index_buf[i] || !s->motion_val_buf[i])
                break;
            s->cur_pic.ref_index[i]  = s->ref_index_buf[i]->data;
            s->cur_pic.motion_val[i] = (int16_t (*)[2])s->motion_val_buf[i]->data + 4;
        }
        if (i < 2) {
            for (i = 0; i < 2; i++) {
                av_buffer_unref(&s->ref_index_buf[i]);
                av_buffer_unref(&s->motion_val_buf[i]);
                s->cur_pic.ref_index[i]  = NULL;
                s->cur_pic.motion_val[i] = NULL;
            }
            return;
        }
    }

    if (s->avctx->debug & FF_DEBUG_ER) {
        for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
            for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
                int status = s->error_status_table[mb_x + mb_y * s->mb_stride];

                av_log(s->avctx, AV_LOG_DEBUG, "%2X ", status);
            }
            av_log(s->avctx, AV_LOG_DEBUG, "\n");
        }
    }

#if 1
    /* handle overlapping slices */
    for (error_type = 1; error_type <= 3; error_type++) {
        int end_ok = 0;

        for (i = s->mb_num - 1; i >= 0; i--) {
            const int mb_xy = s->mb_index2xy[i];
            int error       = s->error_status_table[mb_xy];

            // if error == 2, 4, 8, 16, 32, 64
            if (error & (1 << error_type))
                end_ok = 1;
            if (error & (8 << error_type))
                end_ok = 1;

            if (!end_ok)
                s->error_status_table[mb_xy] |= 1 << error_type;

            if (error & VP_START)
                end_ok = 0;
        }
    }
#endif
#if 1
    /* handle slices with partitions of different length */
    if (s->partitioned_frame) {
        int end_ok = 0;

        for (i = s->mb_num - 1; i >= 0; i--) {
            const int mb_xy = s->mb_index2xy[i];
            int error       = s->error_status_table[mb_xy];

            if (error & ER_AC_END)
                end_ok = 0;
            if ((error & ER_MV_END) ||
                (error & ER_DC_END) ||
                (error & ER_AC_ERROR))
                end_ok = 1;

            if (!end_ok)
                s->error_status_table[mb_xy]|= ER_AC_ERROR;

            if (error & VP_START)
                end_ok = 0;
        }
    }
#endif
    /* handle missing slices */
    if (s->avctx->err_recognition & AV_EF_EXPLODE) {
        int end_ok = 1;

        // FIXME + 100 hack
        for (i = s->mb_num - 2; i >= s->mb_width + 100; i--) {
            const int mb_xy = s->mb_index2xy[i];
            int error1 = s->error_status_table[mb_xy];
            int error2 = s->error_status_table[s->mb_index2xy[i + 1]];

            if (error1 & VP_START)
                end_ok = 1;

            if (error2 == (VP_START | ER_MB_ERROR | ER_MB_END) &&
                error1 != (VP_START | ER_MB_ERROR | ER_MB_END) &&
                ((error1 & ER_AC_END) || (error1 & ER_DC_END) ||
                (error1 & ER_MV_END))) {
                // end & uninit
                end_ok = 0;
            }

            if (!end_ok)
                s->error_status_table[mb_xy] |= ER_MB_ERROR;
        }
    }

#if 1
    /* backward mark errors */
    distance = 9999999;
    for (error_type = 1; error_type <= 3; error_type++) {
        for (i = s->mb_num - 1; i >= 0; i--) {
            const int mb_xy = s->mb_index2xy[i];
            int       error = s->error_status_table[mb_xy];

            if (!s->mbskip_table || !s->mbskip_table[mb_xy]) // FIXME partition specific
                distance++;
            if (error & (1 << error_type))
                distance = 0;

            if (s->partitioned_frame) {
                if (distance < threshold_part[error_type - 1])
                    s->error_status_table[mb_xy] |= 1 << error_type;
            } else {
                if (distance < threshold)
                    s->error_status_table[mb_xy] |= 1 << error_type;
            }

            if (error & VP_START)
                distance = 9999999;
        }
    }
#endif

    /* forward mark errors */
    error = 0;
    for (i = 0; i < s->mb_num; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int old_error   = s->error_status_table[mb_xy];

        if (old_error & VP_START) {
            error = old_error & ER_MB_ERROR;
        } else {
            error |= old_error & ER_MB_ERROR;
            s->error_status_table[mb_xy] |= error;
        }
    }
#if 1
    /* handle not partitioned case */
    if (!s->partitioned_frame) {
        for (i = 0; i < s->mb_num; i++) {
            const int mb_xy = s->mb_index2xy[i];
            int error = s->error_status_table[mb_xy];
            if (error & ER_MB_ERROR)
                error |= ER_MB_ERROR;
            s->error_status_table[mb_xy] = error;
        }
    }
#endif

    dc_error = ac_error = mv_error = 0;
    for (i = 0; i < s->mb_num; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int error = s->error_status_table[mb_xy];
        if (error & ER_DC_ERROR)
            dc_error++;
        if (error & ER_AC_ERROR)
            ac_error++;
        if (error & ER_MV_ERROR)
            mv_error++;
    }
    av_log(s->avctx, AV_LOG_INFO, "concealing %d DC, %d AC, %d MV errors in %c frame\n",
           dc_error, ac_error, mv_error, av_get_picture_type_char(s->cur_pic.f->pict_type));

    s->cur_pic.f->decode_error_flags |= FF_DECODE_ERROR_CONCEALMENT_ACTIVE;

    is_intra_likely = is_intra_more_likely(s);
    /* set unknown mb-type to most likely */
    for (i = 0; i < s->mb_num; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int error = s->error_status_table[mb_xy];
        if (!((error & ER_DC_ERROR) && (error & ER_MV_ERROR)))
            continue;

        if (is_intra_likely)
            s->cur_pic.mb_type[mb_xy] = MB_TYPE_INTRA4x4;
        else
            s->cur_pic.mb_type[mb_xy] = MB_TYPE_16x16 | MB_TYPE_L0;
    }

    // change inter to intra blocks if no reference frames are available
    if (!(s->last_pic.f && s->last_pic.f->data[0]) &&
        !(s->next_pic.f && s->next_pic.f->data[0]))
        for (i = 0; i < s->mb_num; i++) {
            const int mb_xy = s->mb_index2xy[i];
            if (!IS_INTRA(s->cur_pic.mb_type[mb_xy]))
                s->cur_pic.mb_type[mb_xy] = MB_TYPE_INTRA4x4;
        }


    /* handle inter blocks with damaged AC */
    for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
            const int mb_xy   = mb_x + mb_y * s->mb_stride;
            const int mb_type = s->cur_pic.mb_type[mb_xy];
            const int dir     = !(s->last_pic.f && s->last_pic.f->data[0]);
            const int mv_dir  = dir ? MV_DIR_BACKWARD : MV_DIR_FORWARD;
            int mv_type;

            int error = s->error_status_table[mb_xy];

            if (IS_INTRA(mb_type))
                continue; // intra
            if (error & ER_MV_ERROR)
                continue; // inter with damaged MV
            if (!(error & ER_AC_ERROR))
                continue; // undamaged inter

            if (IS_8X8(mb_type)) {
                int mb_index = mb_x * 2 + mb_y * 2 * s->b8_stride;
                int j;
                mv_type = MV_TYPE_8X8;
                for (j = 0; j < 4; j++) {
                    s->mv[0][j][0] = s->cur_pic.motion_val[dir][mb_index + (j & 1) + (j >> 1) * s->b8_stride][0];
                    s->mv[0][j][1] = s->cur_pic.motion_val[dir][mb_index + (j & 1) + (j >> 1) * s->b8_stride][1];
                }
            } else {
                mv_type     = MV_TYPE_16X16;
                s->mv[0][0][0] = s->cur_pic.motion_val[dir][mb_x * 2 + mb_y * 2 * s->b8_stride][0];
                s->mv[0][0][1] = s->cur_pic.motion_val[dir][mb_x * 2 + mb_y * 2 * s->b8_stride][1];
            }

            s->decode_mb(s->opaque, 0 /* FIXME H.264 partitioned slices need this set */,
                         mv_dir, mv_type, &s->mv, mb_x, mb_y, 0, 0);

        }
    }
    // kjlee: run edge detection for full frame 
    edge_detection_picture(s, s->mb_width, s->mb_height, s->mb_stride);

    
    /* guess MVs */
    if (s->cur_pic.f->pict_type == AV_PICTURE_TYPE_B) {
        for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
            for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
                int       xy      = mb_x * 2 + mb_y * 2 * s->b8_stride;
                const int mb_xy   = mb_x + mb_y * s->mb_stride;
                const int mb_type = s->cur_pic.mb_type[mb_xy];
                int mv_dir = MV_DIR_FORWARD | MV_DIR_BACKWARD;

                int error = s->error_status_table[mb_xy];

                if (IS_INTRA(mb_type))
                    continue;
                if (!(error & ER_MV_ERROR))
                    continue; // inter with undamaged MV
                if (!(error & ER_AC_ERROR))
                    continue; // undamaged inter

                if (!(s->last_pic.f && s->last_pic.f->data[0]))
                    mv_dir &= ~MV_DIR_FORWARD;
                if (!(s->next_pic.f && s->next_pic.f->data[0]))
                    mv_dir &= ~MV_DIR_BACKWARD;

                if (s->pp_time) {
                    int time_pp = s->pp_time;
                    int time_pb = s->pb_time;

                    av_assert0(s->avctx->codec_id != AV_CODEC_ID_H264);
                    ff_thread_await_progress(s->next_pic.tf, mb_y, 0);

                    s->mv[0][0][0] = s->next_pic.motion_val[0][xy][0] *  time_pb            / time_pp;
                    s->mv[0][0][1] = s->next_pic.motion_val[0][xy][1] *  time_pb            / time_pp;
                    s->mv[1][0][0] = s->next_pic.motion_val[0][xy][0] * (time_pb - time_pp) / time_pp;
                    s->mv[1][0][1] = s->next_pic.motion_val[0][xy][1] * (time_pb - time_pp) / time_pp;
                } else {
                    s->mv[0][0][0] = 0;
                    s->mv[0][0][1] = 0;
                    s->mv[1][0][0] = 0;
                    s->mv[1][0][1] = 0;
                }

                s->decode_mb(s->opaque, 0, mv_dir, MV_TYPE_16X16, &s->mv,
                             mb_x, mb_y, 0, 0);
            }
        }
    } else{
        //kjlee fix
        //guess_mv(s);
        guess_mv_edge(s);
    }
        

    /* fill DC for inter blocks */
    for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
            int dc, dcu, dcv, y, n;
            int16_t *dc_ptr;
            uint8_t *dest_y, *dest_cb, *dest_cr;
            const int mb_xy   = mb_x + mb_y * s->mb_stride;
            const int mb_type = s->cur_pic.mb_type[mb_xy];

            // error = s->error_status_table[mb_xy];

            if (IS_INTRA(mb_type) && s->partitioned_frame)
                continue;
            // if (error & ER_MV_ERROR)
            //     continue; // inter data damaged FIXME is this good?
            dest_y  = s->cur_pic.f->data[0] + mb_x * 16 + mb_y * 16 * linesize[0];
            dest_cb = s->cur_pic.f->data[1] + mb_x *  8 + mb_y *  8 * linesize[1];
            dest_cr = s->cur_pic.f->data[2] + mb_x *  8 + mb_y *  8 * linesize[2];

            dc_ptr = &s->dc_val[0][mb_x * 2 + mb_y * 2 * s->b8_stride];
            for (n = 0; n < 4; n++) {
                dc = 0;
                for (y = 0; y < 8; y++) {
                    int x;
                    for (x = 0; x < 8; x++)
                       dc += dest_y[x + (n & 1) * 8 +
                             (y + (n >> 1) * 8) * linesize[0]];
                }
                dc_ptr[(n & 1) + (n >> 1) * s->b8_stride] = (dc + 4) >> 3;
            }

            if (!s->cur_pic.f->data[2])
                continue;

            dcu = dcv = 0;
            for (y = 0; y < 8; y++) {
                int x;
                for (x = 0; x < 8; x++) {
                    dcu += dest_cb[x + y * linesize[1]];
                    dcv += dest_cr[x + y * linesize[2]];
                }
            }
            s->dc_val[1][mb_x + mb_y * s->mb_stride] = (dcu + 4) >> 3;
            s->dc_val[2][mb_x + mb_y * s->mb_stride] = (dcv + 4) >> 3;
        }
    }
#if 1
    /* guess DC for damaged blocks */
    guess_dc(s, s->dc_val[0], s->mb_width*2, s->mb_height*2, s->b8_stride, 1);
    guess_dc(s, s->dc_val[1], s->mb_width  , s->mb_height  , s->mb_stride, 0);
    guess_dc(s, s->dc_val[2], s->mb_width  , s->mb_height  , s->mb_stride, 0);

    //kjlee
    //guess_dominant_edge(s, s->mb_width, s->mb_height, s->mb_stride);
    if(temp ==0){
         
         char buf[1024];
         snprintf(buf, sizeof(buf), "%s-%d-lum.ppm", "test", 0);
         FILE *fp = fopen(buf, "w+");
        fprintf(fp, "P6\n%d %d\n255\n", s->cur_pic.f->width, s->cur_pic.f->height); //  The PPM file adds fixed header information.

        uint8_t* temp_edge = (uint8_t*)malloc(s->cur_pic.f->width * s->cur_pic.f->height*3);
        for(int k = 0 ; k < s->cur_pic.f->width * s->cur_pic.f->height ; k++){
            temp_edge[3*k] = (uint8_t)s->cur_pic.edge_mag[k];
            temp_edge[3*k+1] = (uint8_t)s->cur_pic.edge_mag[k];
            temp_edge[3*k+2] = (uint8_t)s->cur_pic.edge_mag[k];
        }

    for (int y=0; y<s->cur_pic.f->height; y++)
        fwrite(&temp_edge[y * s->cur_pic.f->width * 3], 1, s->cur_pic.f->width*3, fp); //ppm storage format
    fclose(fp);
    free(temp_edge);
    temp++;
    }
#endif

    /* filter luma DC */
    filter181(s->dc_val[0], s->mb_width * 2, s->mb_height * 2, s->b8_stride);

#if 1
    /* render DC only intra */
    /*for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
            uint8_t *dest_y, *dest_cb, *dest_cr;
            const int mb_xy   = mb_x + mb_y * s->mb_stride;
            const int mb_type = s->cur_pic.mb_type[mb_xy];

            int error = s->error_status_table[mb_xy];

            if (IS_INTER(mb_type))
                continue;
            if (!(error & ER_AC_ERROR))
                continue; // undamaged

            dest_y  = s->cur_pic.f->data[0] + mb_x * 16 + mb_y * 16 * linesize[0];
            dest_cb = s->cur_pic.f->data[1] + mb_x *  8 + mb_y *  8 * linesize[1];
            dest_cr = s->cur_pic.f->data[2] + mb_x *  8 + mb_y *  8 * linesize[2];
            if (!s->cur_pic.f->data[2])
                dest_cb = dest_cr = NULL;

            put_dc(s, dest_y, dest_cb, dest_cr, mb_x, mb_y);
            
        }
    }*/

    //kjlee
     uint8_t *dest_y, *dest_cb, *dest_cr;
    dest_y  = s->cur_pic.f->data[0];
    dest_cb = s->cur_pic.f->data[1];
    dest_cr = s->cur_pic.f->data[2];
    for (mb_y = 0; mb_y < s->mb_height; mb_y++) {
        for (mb_x = 0; mb_x < s->mb_width; mb_x++) {
           
            const int mb_xy   = mb_x + mb_y * s->mb_stride;
            const int mb_type = s->cur_pic.mb_type[mb_xy];

            int error = s->error_status_table[mb_xy];

            if (IS_INTER(mb_type))
                continue;
            if (!(error & ER_AC_ERROR))
                continue; // undamaged

            
            if (!s->cur_pic.f->data[2])
                dest_cb = dest_cr = NULL;

            //put_dc(s, dest_y, dest_cb, dest_cr, mb_x, mb_y);
            //kjlee
            put_dc_edge(s, dest_y, dest_cb, dest_cr, mb_x, mb_y, s->mb_stride);

        }
    }
#endif

    if (s->avctx->error_concealment & FF_EC_DEBLOCK) {
        /* filter horizontal block boundaries */
        h_block_filter(s, s->cur_pic.f->data[0], s->mb_width * 2,
                       s->mb_height * 2, linesize[0], 1);

        /* filter vertical block boundaries */
        v_block_filter(s, s->cur_pic.f->data[0], s->mb_width * 2,
                       s->mb_height * 2, linesize[0], 1);

        if (s->cur_pic.f->data[2]) {
            h_block_filter(s, s->cur_pic.f->data[1], s->mb_width,
                        s->mb_height, linesize[1], 0);
            h_block_filter(s, s->cur_pic.f->data[2], s->mb_width,
                        s->mb_height, linesize[2], 0);
            v_block_filter(s, s->cur_pic.f->data[1], s->mb_width,
                        s->mb_height, linesize[1], 0);
            v_block_filter(s, s->cur_pic.f->data[2], s->mb_width,
                        s->mb_height, linesize[2], 0);
        }
    }

    /* clean a few tables */
    for (i = 0; i < s->mb_num; i++) {
        const int mb_xy = s->mb_index2xy[i];
        int       error = s->error_status_table[mb_xy];

        if (s->mbskip_table && s->cur_pic.f->pict_type != AV_PICTURE_TYPE_B &&
            (error & (ER_DC_ERROR | ER_MV_ERROR | ER_AC_ERROR))) {
            s->mbskip_table[mb_xy] = 0;
        }
        if (s->mbintra_table)
            s->mbintra_table[mb_xy] = 1;
    }

    for (i = 0; i < 2; i++) {
        av_buffer_unref(&s->ref_index_buf[i]);
        av_buffer_unref(&s->motion_val_buf[i]);
        s->cur_pic.ref_index[i]  = NULL;
        s->cur_pic.motion_val[i] = NULL;
    }

    memset(&s->cur_pic, 0, sizeof(ERPicture));
    memset(&s->last_pic, 0, sizeof(ERPicture));
    memset(&s->next_pic, 0, sizeof(ERPicture));
}
