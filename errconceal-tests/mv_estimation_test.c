/*
 * Copyright (c) 2001 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * @file
 * video decoding with libavcodec API example
 *
 * @example decode_video.c
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libavcodec/avcodec.h>
#include <libswscale/swscale.h>
#include <libavutil/imgutils.h>
#include <libavutil/motion_vector.h>

#define INBUF_SIZE 4096


static struct SwsContext *sws_ctx = NULL;
typedef struct Decoder{
    const AVCodec *codec;
    AVCodecParserContext *parser;
    AVCodecContext *context ;
    AVFrame *frame;
    AVFrame *frame_rgb;
    AVPacket *pkt;
    uint8_t inbuf[INBUF_SIZE + AV_INPUT_BUFFER_PADDING_SIZE];
    uint8_t *data;
    size_t data_size;
    AVDictionary *opts;
    int ret;
    int *mv;
    int mv_size;

} Decoder;

static void pgm_save(unsigned char *buf, int wrap, int xsize, int ysize,
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
void ppm_save_raw(uint8_t* data, int width, int height, char *filename){
    int y;
    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "P6\n%d %d\n255\n", width, height); //  The PPM file adds fixed header information.

    for (y=0; y<height; y++)
        fwrite(data+y*width*3, 1, width*3, fp); //ppm storage format
    fclose(fp);
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

static void decode(AVCodecContext *dec_ctx, AVFrame *frame, AVFrame *frame_rgb, AVPacket *pkt, Decoder* dec)
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
        //snprintf(buf, sizeof(buf), "%s-%d.ppm", filename, dec_ctx->frame_number);
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
        frame_rgb->pict_type = frame->pict_type;
        sws_scale(sws_ctx, frame->data, frame->linesize,
                 0, dec_ctx->height, frame_rgb->data, frame_rgb->linesize);
    
        if (frame->pict_type == AV_PICTURE_TYPE_P) {
            AVFrameSideData *sd;
            sd = av_frame_get_side_data(frame, AV_FRAME_DATA_MOTION_VECTORS);
            if (sd) {
              const AVMotionVector *mvs = (const AVMotionVector *)sd->data;
              int num_data = sd->size / sizeof(*mvs);
              dec->mv_size = num_data;
              dec->mv = (int*) malloc(sizeof(int) * num_data * 2);
              for (int i = 0; i < sd->size / sizeof(*mvs); i++) {
                    const AVMotionVector *mv = &mvs[i];
                    dec->mv[2*i] = mv->motion_x;
                    dec->mv[2*i+1] = mv->motion_y;
                }
            }
        }   
    

    
        //    ppm_save(frame_rgb, frame_rgb->width, frame_rgb->height, buf);

    }
}

void initDecoder(Decoder *dec)
{
    dec->pkt = av_packet_alloc();
    if (!dec->pkt)
        exit(1);
    dec->context = NULL;
    dec->opts = NULL;
    // set end of buffer to 0 (this ensures that no overreading happens for damaged MPEG streams) 
    memset(dec->inbuf + INBUF_SIZE, 0, AV_INPUT_BUFFER_PADDING_SIZE);

    // find the h264 video decoder 
    dec->codec = avcodec_find_decoder(AV_CODEC_ID_H264);
    if (!dec->codec) {
        fprintf(stderr, "Codec not found\n");
        exit(1);
    }

    dec->parser = av_parser_init(dec->codec->id);
    if (!dec->parser) {
        fprintf(stderr, "parser not found\n");
        exit(1);
    }

    dec->context = avcodec_alloc_context3(dec->codec);
    if (!dec->context) {
        fprintf(stderr, "Could not allocate video codec context\n");
        exit(1);
    }

    av_dict_set(&dec->opts, "flags2", "+export_mvs", 0);

    if (avcodec_open2(dec->context, dec->codec, &dec->opts) < 0) {
        fprintf(stderr, "Could not open codec\n");
        exit(1);
    }


    dec->frame = av_frame_alloc();
    if (!dec->frame) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }

    dec->frame_rgb = av_frame_alloc();
    if (!dec->frame_rgb) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }



}



void decodeFrame(Decoder *dec, FILE *f)
{

    char isFrame = 0;
    while(1){
        
        dec->data_size = fread(dec->inbuf, 1, INBUF_SIZE, f);
        if (!dec->data_size) break;
        dec->data = dec->inbuf;
        while (dec->data_size > 0) {


            dec->ret = av_parser_parse2(dec->parser, dec->context, &dec->pkt->data, &dec->pkt->size,
                                    dec->data, dec->data_size, AV_NOPTS_VALUE, AV_NOPTS_VALUE, 0);

            if (dec->ret < 0) {
                fprintf(stderr, "Error while parsing\n");
                exit(1);
            }
            dec->data      += dec->ret;
            dec->data_size -= dec->ret;
            if (dec->pkt->size){
                decode(dec->context, dec->frame, dec->frame_rgb, dec->pkt, dec);
                isFrame = 1;
            }
                
        }
        if(isFrame) break;

    }
    
}


void freeDecoder(Decoder *dec)
{
    av_parser_close(dec->parser);
    avcodec_free_context(&dec->context);
    av_frame_free(&dec->frame);
    av_frame_free(&dec->frame_rgb);
    av_packet_free(&dec->pkt);

}

void compareMV(Decoder *gt_dec, Decoder *loss_dec, char* filename)
{    

    // score using motion vectors
    
    if (gt_dec->frame_rgb->pict_type != AV_PICTURE_TYPE_P) {
        printf("gt frame not ptype\n");
        return;
    }
    if (loss_dec->frame_rgb->pict_type != AV_PICTURE_TYPE_P) return;
   


    uint8_t* output = (uint8_t*) malloc(gt_dec->mv_size);
    memset(output, 0, gt_dec->mv_size);

    for (int i = 0; i < gt_dec->mv_size ; i++) {
        float mv_x = (float) gt_dec->mv[2*i];
        float mv_y = (float) gt_dec->mv[2*i+1];

        float pred_mv_x = (float) loss_dec->mv[2*i];
        float pred_mv_y = (float) loss_dec->mv[2*i+1];

        
        float dist = (mv_x - pred_mv_x)*(mv_x - pred_mv_x) + (mv_y - pred_mv_y) * (mv_y - pred_mv_y);
        dist = sqrt(dist);
        if(dist > 1000) output[i] = 255;
    }

    int tot_mb_width = gt_dec->frame_rgb->width / 16;
    printf("tot_mb_width : %d %d\n", tot_mb_width, gt_dec->mv_size);

    uint8_t* score = (uint8_t*) malloc(gt_dec->frame_rgb->height * gt_dec->frame_rgb->width * 3);
    for(int i = 0 ; i < gt_dec->frame_rgb->height ; i++){
        for(int j = 0 ; j < gt_dec->frame_rgb->width ; j++){
            int mb_height = i / 16;
            int mb_width = j / 16;
            if(output[mb_height * tot_mb_width + mb_width] == 255) printf("test\n");
            score[3*(i * gt_dec->frame_rgb->width + j)] = output[mb_height * tot_mb_width + mb_width]; 
            score[3*(i * gt_dec->frame_rgb->width + j)+1] = output[mb_height * tot_mb_width + mb_width]; 
            score[3*(i * gt_dec->frame_rgb->width + j)+2] = output[mb_height * tot_mb_width + mb_width]; 
        }
    }

    ppm_save_raw(score, gt_dec->frame_rgb->width, gt_dec->frame_rgb->height, filename);




  
    
}

int main(int argc, char **argv)
{
    const char *filename, *lossfilename, *outfilename, *scorefilename;
    FILE *f, *lossf;

    Decoder* dec = (Decoder*) malloc(sizeof(Decoder));
    Decoder* lossdec = (Decoder*) malloc(sizeof(Decoder));

    
    if (argc <= 4) {
        fprintf(stderr, "Usage: %s <input file> <output file>\n"
                "And check your input file is encoded by mpeg1video please.\n", argv[0]);
        exit(0);
    }
    filename    = argv[1];
    lossfilename = argv[2];
    outfilename = argv[3];
    scorefilename = argv[4];


    f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }

    lossf = fopen(lossfilename, "rb");
    if(!lossf){
        fprintf(stderr, "Could not open %s\n", lossfilename);
        exit(1);
    }

    initDecoder(dec);
    initDecoder(lossdec);
    printf("init Decoder\n");


    int cnt = 0;    
    char buf[1024];

    while (1) {
        if(feof(f) || feof(lossf)) break;
        
        decodeFrame(dec, f);
        decodeFrame(lossdec, lossf);
        snprintf(buf, sizeof(buf), "%s-%d.ppm", scorefilename, cnt);

        compareMV(dec, lossdec, buf);
        cnt++;
        if(cnt == 10) break;
    }

    /* flush the decoder */
    decode(dec->context, dec->frame, dec->frame_rgb, NULL, dec);
    decode(lossdec->context, lossdec->frame, lossdec->frame_rgb, NULL, dec);


    fclose(f);
    fclose(lossf);
    freeDecoder(dec);
    freeDecoder(lossdec);
    free(dec);
    free(lossdec);

    return 0;
}
