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

typedef struct {
	int width;
	int height;
	int **imageData;
	int **gx;
	int **gy;
} Edge;
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


void sobel_edge_detector(uint8_t* data, Edge* output)
{
    int i, j, gx, gy;
	int mx[3][3] = {
		{-1, 0, 1},
		{-2, 0, 2},
		{-1, 0, 1}
	};
	int my[3][3] = {
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

void edge_save(uint8_t* data, int width, int height, char *filename)
{
    Edge out;
    out.width = width;
    out.height = height;
    out.imageData = (int**) calloc(out.height, sizeof(int*));
    int i,j;
	for(i = 0; i < out.height; i++) {
		out.imageData[i] = calloc(out.width, sizeof(int));
	}
	
	out.gx = (int**) calloc(out.height, sizeof(int*));
	for(i = 0; i < out.height; i++) {
		out.gx[i] = calloc(out.width, sizeof(int));
	}
	
	out.gy = (int**) calloc(out.height, sizeof(int*));
	for(i = 0; i < out.height; i++) {
		out.gy[i] = calloc(out.width, sizeof(int));
	}

    for(i = 0; i < out.height; i++) {
		for(j = 0; j < out.width; j++) {
			out.imageData[i][j] = data[i*width+j];
			out.gx[i][j] = data[i*width+j];
			out.gy[i][j] = data[i*width+j];
		};
	}
    sobel_edge_detector(data, &out);
    uint8_t* edge_data = (uint8_t*) malloc(sizeof(uint8_t)*width*height);
    for(i = 0 ; i < height ; i++){
        for(j = 0 ; j < width ; j++){
            edge_data[i*width+j] = out.imageData[i][j];
        }
    }
    pgm_save(edge_data, width, width, height, filename);

    for(int i = 0 ; i < out.height; i++){
            free(out.imageData[i]);
            free(out.gx[i]);
            free(out.gy[i]);
    }
    free(out.imageData);
    free(out.gx);
    free(out.gy);

   
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
        sws_scale(sws_ctx, frame->data, frame->linesize,
                 0, dec_ctx->height, frame_rgb->data, frame_rgb->linesize);


        snprintf(buf, sizeof(buf), "%s-%d-lum.ppm", filename, dec_ctx->frame_number);
        uint8_t* luminance = (uint8_t*) malloc (sizeof(uint8_t) * frame->width * frame->height);
        yuvframe_to_luminance(frame, luminance);
        edge_save(luminance, frame->width, frame->height, buf);

        snprintf(buf, sizeof(buf), "%s-%d-grey.ppm", filename, dec_ctx->frame_number);
        uint8_t* grey = (uint8_t*) malloc(sizeof(uint8_t) * frame_rgb->width * frame_rgb->height);
        rgbframe_to_grey(frame_rgb, grey);
        edge_save(grey, frame_rgb->width, frame_rgb->height, buf);
        free(luminance);
        free(grey);

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
        if(cnt == 30) break;
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