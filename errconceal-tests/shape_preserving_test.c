
/**
 * @file
 * sobel edge detector test code
 *
 * @example edge_detector_test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libavcodec/avcodec.h>
#include <libswscale/swscale.h>
#include <libavutil/imgutils.h>

#define INBUF_SIZE 4096
#define MB_SIZE 16
#define MB_ERROR_WIDTH_START 49
#define MB_ERROR_WIDTH_END 50
#define MB_ERROR_HEIGHT_START 8
#define MB_ERROR_HEIGHT_END 10

#define LEFT 0
#define RIGHT 1
#define TOP 2
#define BOTTOM 3


typedef struct {
	int width; // with padding for macroblock
	int height;
    int orig_width;
    int orig_height;
    int mb_width;
    int mb_height;
    float *magnitude;
    float *direction;
    int **mb_error;
	int **imageData;
	int **gx;
	int **gy;
} Edge;

typedef struct{
    float magnitude;
    float direction;

} EdgePixel;

static struct SwsContext *sws_ctx = NULL;

static void pgm_save(uint8_t *buf, int wrap, int xsize, int ysize,
                     char *filename)
{
    FILE *f;
    int i;

    f = fopen(filename,"wb");
    fprintf(f, "P5\n%d %d\n%d\n", xsize, ysize, 255);
    for (i = 0; i < ysize; i++)
        fwrite(buf + i * wrap, 1, xsize, f);

    fclose(f);
}

void ppm_save(AVFrame *pframe, int width, int height, char *filename)
{
    int y;
    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "P6\n%d %d\n255\n", width, height); //  The PPM file adds fixed header information.

    for (y=0; y<height; y++)
        fwrite(pframe->data[0]+y*pframe->linesize[0], 1, width*3, fp); //ppm storage format
    fclose(fp);
}

void yuvframe_to_luminance(AVFrame* yuv, uint8_t* out){
    for(int i = 0 ; i < yuv->height ; i++){
        for(int j = 0 ; j < yuv->width ; j++){
            out[i*yuv->width + j] = yuv->data[0][i*yuv->linesize[0] + j];
        }
    }
}

void conceal_yuv_with_edge(AVFrame* yuv, Edge* data)
{
    float x_inc = 0;
    float y_inc = 0;
    int *reference_pixels = calloc(yuv->width * yuv->height * 4, sizeof(int));

    for(int i = 0 ; i < yuv->height ; i++){
        for(int j = 0 ; j < yuv->width  ;j++){
            int mb_height_idx = i / MB_SIZE;
            int mb_width_idx = j / MB_SIZE;
            if(data->mb_error[mb_height_idx][mb_width_idx] == 0 ) continue;
            float direction = tan(data->direction[mb_height_idx * data->mb_width + mb_width_idx]);

           // printf("theta: %f\n",data->direction[mb_height_idx * data->mb_width + mb_width_idx]*180/M_PI);
           // printf("direction : %f\n", direction);
           // printf("mag: %f\n", data->magnitude[mb_height_idx * data->mb_width + mb_width_idx]);

            if(abs(direction) >= 1){
                x_inc = 1 / direction;
                y_inc = 1;
            } else{
                x_inc = 1;
                y_inc = direction;
            }

            float curr_x = (float) j;
            float curr_y = (float) i;
            int curr_pix_x = round(curr_x);
            int curr_pix_y = round(curr_y);
            int curr_mb_x = curr_pix_x / MB_SIZE;
            int curr_mb_y = curr_pix_y / MB_SIZE;
            while(1){
                curr_x += x_inc;
                curr_y -= y_inc;

                curr_pix_x = round(curr_x);
                curr_pix_y = round(curr_y);

                curr_mb_x = curr_pix_x / MB_SIZE;
                curr_mb_y = curr_pix_y / MB_SIZE;
                if(data->mb_error[curr_mb_y][curr_mb_x] == 0){
                    reference_pixels[4*(i*yuv->width + j)+0] = curr_pix_x;
                    reference_pixels[4*(i*yuv->width + j)+1] = curr_pix_y;

                    break;
                }

                if(curr_pix_x >= yuv->width || curr_pix_x < 0 || curr_pix_y >= yuv->height || curr_pix_y < 0) break;
            }
            while(1){
                curr_x -= x_inc;
                curr_y += y_inc;

                curr_pix_x = round(curr_x);
                curr_pix_y = round(curr_y);

                curr_mb_x = curr_pix_x / MB_SIZE;
                curr_mb_y = curr_pix_y / MB_SIZE;
                if(data->mb_error[curr_mb_y][curr_mb_x] == 0){
                    reference_pixels[4*(i*yuv->width + j)+2] = curr_pix_x;
                    reference_pixels[4*(i*yuv->width + j)+3] = curr_pix_y;
                    break;
                }

                if(curr_pix_x >= yuv->width || curr_pix_x < 0 || curr_pix_y >= yuv->height || curr_pix_y < 0) break;
            }
        }
    }

    for(int i = 0 ; i < yuv->height ; i++){
        for(int j = 0 ; j < yuv->width  ;j++){
            int mb_height_idx = i / MB_SIZE;
            int mb_width_idx = j / MB_SIZE;
            if(data->mb_error[mb_height_idx][mb_width_idx] == 0 ) continue;

            if(reference_pixels[4*(i*yuv->width + j)+0] < 0 || reference_pixels[4*(i*yuv->width + j)+3] < 0){
                yuv->data[0][i*yuv->linesize[0] + j] = 0;
                yuv->data[1][i/2*yuv->linesize[1] + j/2] = 0;
                yuv->data[2][i/2*yuv->linesize[2] + j/2] = 0;
                continue;
            }

            int right_x = reference_pixels[4*(i*yuv->width + j)+0];
            int right_y = reference_pixels[4*(i*yuv->width + j)+1];
            
            int left_x = reference_pixels[4*(i*yuv->width + j)+2];
            int left_y = reference_pixels[4*(i*yuv->width + j)+3];
            float dist_right = sqrt((j-right_x)*(j-right_x)+(i-right_y)*(i-right_y));
            float dist_left = sqrt((j-left_x)*(j-left_x) + (i-left_y)*(i-left_y));

            yuv->data[0][i*yuv->linesize[0] + j] = round((dist_right * yuv->data[0][left_y*yuv->linesize[0] + left_x] + dist_left * yuv->data[0][right_y*yuv->linesize[0] + right_x]) / (dist_right + dist_left));
            yuv->data[1][i/2*yuv->linesize[1] + j/2] = round((dist_right * yuv->data[1][left_y/2*yuv->linesize[1] + left_x/2] + dist_left * yuv->data[1][right_y/2*yuv->linesize[1] + right_x/2]) / (dist_right + dist_left));
            yuv->data[2][i/2*yuv->linesize[2] + j/2] = round((dist_right * yuv->data[2][left_y/2*yuv->linesize[2] + left_x/2] + dist_left * yuv->data[2][right_y/2*yuv->linesize[2] + right_x/2]) / (dist_right + dist_left));
           // todo: move this to save both
           // this is for testing only
           // yuv->data[0][i*yuv->linesize[0] + j] = 0;
           // yuv->data[1][i/2*yuv->linesize[1] + j/2] = 0;
           // yuv->data[2][i/2*yuv->linesize[2] + j/2] = 0;
        
           
        }
    
    }
    free(reference_pixels);
}

void rgbframe_to_grey(AVFrame* rgb, uint8_t* out)
{

    for(int y = 0 ; y < rgb->height ; y++){
        for(int x = 0 ; x < rgb->width ; x++){
            int curr = y*rgb->linesize[0] + x*3;
            out[y*rgb->width + x] = 0.3*rgb->data[0][curr+0] + 0.59*rgb->data[0][curr+1] + 0.11*rgb->data[0][curr+2];
        }
    }

}

int convolution(Edge* image, int kernel[3][3], int row, int col) {
	int i, j, sum = 0;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			sum += image->imageData[i + row][j + col] * kernel[i][j];
		}
	}
	return sum;
}


void sobel_edge_detector(Edge* output)
{
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
	
	for (i = 1; i < output->height - 2; i++) {
		for (j = 1; j < output->width - 2; j++) {
			gx = convolution(output, mx, i, j);
			gy = convolution(output, my, i, j);
			output->imageData[i][j] = sqrt(gx*gx + gy*gy);
			output->gx[i][j] = gx;
			output->gy[i][j] = gy;
		}
	}
	
}


// add padding to fit macroblocks
void calc_padding(Edge* out)
{
    int remainder = out->orig_width % MB_SIZE;
    int padding = 0;
    if(remainder > 0){
        padding = (out->orig_width / MB_SIZE + 1) * MB_SIZE - out->orig_width;
    }
    out->width = out->orig_width + padding;

    remainder = out->orig_height % MB_SIZE;
    padding = 0;
    if(remainder > 0){
        padding = (out->orig_height / MB_SIZE + 1) * MB_SIZE - out->orig_height;
    }
    out->height = out->orig_height + padding;

    out->mb_width = out->width / MB_SIZE;
    out->mb_height = out->height / MB_SIZE;

}


// get edge direction with maximum magnitude 
// for each lossy mb, get closest macroblock with no error from 4 direction (left, right, bottom, top) 
// from the closest macroblock get the biggest edge direction from the boundary pixels
// example) left macroblock, look at the pixels at the right side boundary, get the largest edge direction
EdgePixel get_dominant_edge(Edge* data, int mb_width_idx, int mb_height_idx, int direction)
{
    float max_magnitude = 0;
    float max_direction = 0;
    if(direction == LEFT){
        for(int i = 0 ; i < MB_SIZE ; i++){
            if(data->imageData[mb_height_idx*MB_SIZE + i][mb_width_idx*MB_SIZE] > max_magnitude){
                max_magnitude = data->imageData[mb_height_idx*MB_SIZE + i][mb_width_idx*MB_SIZE];
                max_direction = atan2(data->gy[mb_height_idx*MB_SIZE + i][mb_width_idx*MB_SIZE], data->gx[mb_height_idx*MB_SIZE + i][mb_width_idx*MB_SIZE]);
            }
        }
    } else if(direction == RIGHT){
        for(int i = 0 ; i < MB_SIZE ; i++){
            if(data->imageData[mb_height_idx*MB_SIZE + i][(mb_width_idx+1)*MB_SIZE-1] > max_magnitude){
                max_magnitude = data->imageData[mb_height_idx*MB_SIZE + i][(mb_width_idx+1)*MB_SIZE-1];
                max_direction = atan2(data->gy[mb_height_idx*MB_SIZE + i][(mb_width_idx+1)*MB_SIZE-1], data->gx[mb_height_idx*MB_SIZE + i][(mb_width_idx+1)*MB_SIZE-1]);
            }
        }
    } else if(direction == TOP){
        for(int i = 0 ; i < MB_SIZE ; i++){
            if(data->imageData[mb_height_idx*MB_SIZE][mb_width_idx*MB_SIZE + i] > max_magnitude){
                max_magnitude = data->imageData[mb_height_idx*MB_SIZE][mb_width_idx*MB_SIZE + i];
                max_direction = atan2(data->gy[mb_height_idx*MB_SIZE][mb_width_idx*MB_SIZE + i], data->gx[mb_height_idx*MB_SIZE][mb_width_idx*MB_SIZE + i]);
            }
        }

    } else if(direction == BOTTOM){
        for(int i = 0 ; i < MB_SIZE ; i++){
            if(data->imageData[(mb_height_idx+1)*MB_SIZE-1][mb_width_idx*MB_SIZE + i] > max_magnitude){
                max_magnitude = data->imageData[(mb_height_idx+1)*MB_SIZE-1][mb_width_idx*MB_SIZE + i];
                max_direction = atan2(data->gy[(mb_height_idx+1)*MB_SIZE-1][mb_width_idx*MB_SIZE + i], data->gx[(mb_height_idx+1)*MB_SIZE-1][mb_width_idx*MB_SIZE + i]);
            }
        }
    }
    EdgePixel pixel ;
    pixel.magnitude = max_magnitude;
    pixel.direction = max_direction;
    return pixel;
}


// get dominant edge direction for each macroblock 
void calculate_mb_edge(Edge* data){
    int mb_width = data->mb_width;
    int mb_height = data->mb_height;
    
    data->magnitude = calloc(mb_width * mb_height, sizeof(float));
    data->direction = calloc(mb_width * mb_height, sizeof(float));
    EdgePixel  pixel;
    float curr_mag, curr_dir;
    for(int i = 0 ; i < mb_height ; i++){
        curr_mag = 0;
        curr_dir = 0;
        for(int j = 0 ; j < mb_width ; j++){
            if(data->mb_error[i][j] == 0){
                pixel = get_dominant_edge(data, j, i, RIGHT);
                curr_mag = pixel.magnitude;
                curr_dir = pixel.direction;
            }

            if(curr_mag > data->magnitude[i*mb_width + j]){
                data->magnitude[i*mb_width + j] = curr_mag;
                data->direction[i*mb_width + j] = curr_dir;
            }

        }
    }
    for(int i = 0 ; i < mb_height ; i++){
        curr_mag = 0;
        curr_dir = 0;
        for(int j = mb_width-1 ; j >=0 ; j--){
            if(data->mb_error[i][j] == 0){
                pixel = get_dominant_edge(data, j, i, LEFT);
                curr_mag = pixel.magnitude;
                curr_dir = pixel.direction;
            }

            if(curr_mag > data->magnitude[i*mb_width + j]){
                data->magnitude[i*mb_width + j] = curr_mag;
                data->direction[i*mb_width + j] = curr_dir;
            }

        }
    }
    for(int j = 0 ; j < mb_width ; j++){
        curr_mag = 0;
        curr_dir = 0;
        for(int i = 0 ; i < mb_height ; i++){
            if(data->mb_error[i][j] == 0){

                pixel = get_dominant_edge(data, j, i, BOTTOM);

                curr_mag = pixel.magnitude;
                curr_dir = pixel.direction;
            }

            if(curr_mag > data->magnitude[i*mb_width + j]){
                data->magnitude[i*mb_width + j] = curr_mag;
                data->direction[i*mb_width + j] = curr_dir;
            }

        }
    }
    for(int j = 0 ; j < mb_width ; j++){
        curr_mag = 0;
        curr_dir = 0;
        for(int i = mb_height - 1 ; i >=0 ; i--){
            if(data->mb_error[i][j] == 0){
                pixel = get_dominant_edge(data, j, i, TOP);
                curr_mag = pixel.magnitude;
                curr_dir = pixel.direction;
            }

            if(curr_mag > data->magnitude[i*mb_width + j]){
                data->magnitude[i*mb_width + j] = curr_mag;
                data->direction[i*mb_width + j] = curr_dir;
            }

        }
    }

}


// calculate edge for each pixel, magnitude, gx, gy
void calculate_edge(uint8_t* data, int width, int height, Edge* out)
{
    out->orig_width = width;
    out->orig_height = height;
    calc_padding(out);
    int i,j;

    out->imageData = (int**) calloc(out->height, sizeof(int*));
    for(i = 0; i < out->height; i++) {
		out->imageData[i] = calloc(out->width, sizeof(int));
	}
	
	out->gx = (int**) calloc(out->height, sizeof(int*));
	for(i = 0; i < out->height; i++) {
		out->gx[i] = calloc(out->width, sizeof(int));
	}

	out->gy = (int**) calloc(out->height, sizeof(int*));
	for(i = 0; i < out->height; i++) {
		out->gy[i] = calloc(out->width, sizeof(int));
	}

    out->mb_error = (int**) calloc(out->mb_height, sizeof(int*));
    for(i = 0 ; i < out->mb_height ; i++){
        out->mb_error[i] = calloc(out->mb_width , sizeof(int));
        memset(out->mb_error[i], 0, out->mb_width * sizeof(int));
    }

    for(i = MB_ERROR_HEIGHT_START ; i < MB_ERROR_HEIGHT_END ; i++){
        for(j = MB_ERROR_WIDTH_START ; j < MB_ERROR_WIDTH_END ; j++){
            out->mb_error[i][j] = 1;
        }
    }

    for(i = 0; i < height; i++) {
		for(j = 0; j < width; j++) {
			out->imageData[i][j] = data[i*width+j];
			out->gx[i][j] = data[i*width+j];
			out->gy[i][j] = data[i*width+j];
		};
	}

    sobel_edge_detector(out);


}


static void decode(AVCodecContext *dec_ctx, AVFrame *frame, AVFrame *frame_rgb, AVPacket *pkt,
                   const char *filename)
{
    char buf[1024];
    int ret;

    ret = avcodec_send_packet(dec_ctx, pkt);
    if (ret < 0) {
        fprintf(stderr, "Error sending a packet for decoding\n");
        exit(1);
    }

    while (ret >= 0) {
        ret = avcodec_receive_frame(dec_ctx, frame);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            return;
        else if (ret < 0) {
            fprintf(stderr, "Error during decoding\n");
            exit(1);
        }
        

        printf("saving frame %3d\n", dec_ctx->frame_number);
        fflush(stdout);

        /* the picture is allocated by the decoder. no need to
           free it */
        //pgm_save(frame->data[0], frame->linesize[0],
         //        frame->width, frame->height, buf);

        if(sws_ctx==NULL){
                sws_ctx = sws_getContext(
                        frame->width,
                        frame->height,
                        AV_PIX_FMT_YUV420P,
                        frame->width,
                        frame->height,
                        AV_PIX_FMT_RGB24,
                        SWS_BILINEAR,
                        NULL,
                        NULL,
                        NULL
                );
                int numBytes = av_image_get_buffer_size(AV_PIX_FMT_RGB24, frame->width, frame->height, 32);
                uint8_t *buffer = (uint8_t *) av_malloc(numBytes * sizeof(uint8_t));
                // TODO: free arrays?
                av_image_fill_arrays(frame_rgb->data, frame_rgb->linesize, buffer, AV_PIX_FMT_RGB24, frame->width, frame->height, 32);
                frame_rgb->width = frame->width;
                frame_rgb->height = frame->height;
                frame_rgb->format = AV_PIX_FMT_RGB24;
        }
       



        snprintf(buf, sizeof(buf), "%s-%d-lum.ppm", filename, dec_ctx->frame_number);
        uint8_t* luminance = calloc (frame->width * frame->height, sizeof(uint8_t) );
        yuvframe_to_luminance(frame, luminance);


        Edge* e = (Edge*) malloc(sizeof(Edge));
        calculate_edge(luminance, frame->width, frame->height, e);
        calculate_mb_edge(e);
        conceal_yuv_with_edge(frame, e);
        
        
        sws_scale(sws_ctx, frame->data, frame->linesize,
                 0, dec_ctx->height, frame_rgb->data, frame_rgb->linesize);
        ppm_save(frame_rgb, frame_rgb->width, frame_rgb->height, buf);

        
        free(luminance);
        free(e->magnitude);
        free(e->direction);
        for(int i = 0 ; i < e->height; i++){
            free(e->imageData[i]);
            free(e->gx[i]);
            free(e->gy[i]);
        }


        for(int i= 0 ; i < e->mb_height; i++){
            free(e->mb_error[i]);
        }
        free(e->imageData);
        free(e->gx);
        free(e->gy);
        free(e->mb_error);

    }

}

int main(int argc, char **argv)
{
    const char *filename, *outfilename;
    const AVCodec *codec;
    AVCodecParserContext *parser;
    AVCodecContext *c= NULL;
    FILE *f;
    AVFrame *frame, *frame_rgb;
    uint8_t inbuf[INBUF_SIZE + AV_INPUT_BUFFER_PADDING_SIZE];
    uint8_t *data;
    size_t   data_size;
    int ret;
    AVPacket *pkt;

    if (argc <= 2) {
        fprintf(stderr, "Usage: %s <input file> <output file>\n"
                "And check your input file is encoded by mpeg1video please.\n", argv[0]);
        exit(0);
    }
    filename    = argv[1];
    outfilename = argv[2];

    pkt = av_packet_alloc();
    if (!pkt)
        exit(1);

    /* set end of buffer to 0 (this ensures that no overreading happens for damaged MPEG streams) */
    memset(inbuf + INBUF_SIZE, 0, AV_INPUT_BUFFER_PADDING_SIZE);

    /* find the h264 video decoder */
    codec = avcodec_find_decoder(AV_CODEC_ID_H264);
    if (!codec) {
        fprintf(stderr, "Codec not found\n");
        exit(1);
    }

    parser = av_parser_init(codec->id);
    if (!parser) {
        fprintf(stderr, "parser not found\n");
        exit(1);
    }

    c = avcodec_alloc_context3(codec);
    if (!c) {
        fprintf(stderr, "Could not allocate video codec context\n");
        exit(1);
    }

    /* For some codecs, such as msmpeg4 and mpeg4, width and height
       MUST be initialized there because this information is not
       available in the bitstream. */

    /* open it */
    if (avcodec_open2(c, codec, NULL) < 0) {
        fprintf(stderr, "Could not open codec\n");
        exit(1);
    }

    f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }

    frame = av_frame_alloc();
    if (!frame) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }

    frame_rgb = av_frame_alloc();
    if (!frame_rgb) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }



    int cnt = 0;
    while (!feof(f)) {
        /* read raw data from the input file */
        data_size = fread(inbuf, 1, INBUF_SIZE, f);
        if (!data_size)
            break;

        /* use the parser to split the data into frames */
        data = inbuf;
        while (data_size > 0) {
            ret = av_parser_parse2(parser, c, &pkt->data, &pkt->size,
                                   data, data_size, AV_NOPTS_VALUE, AV_NOPTS_VALUE, 0);
            if (ret < 0) {
                fprintf(stderr, "Error while parsing\n");
                exit(1);
            }
            data      += ret;
            data_size -= ret;

            if (pkt->size){
                decode(c, frame, frame_rgb, pkt, outfilename);
                cnt++;
            }
                
        }
        if(cnt == 10) break;
    }

    /* flush the decoder */
    decode(c, frame, frame_rgb, NULL, outfilename);

    fclose(f);

    av_parser_close(parser);
    avcodec_free_context(&c);
    av_frame_free(&frame);
    av_frame_free(&frame_rgb);
    av_packet_free(&pkt);

    return 0;
}
