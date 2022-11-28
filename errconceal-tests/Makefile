# use pkg-config for getting CFLAGS and LDLIBS
FFMPEG_LIBS=    libavdevice                        \
                libavformat                        \
                libavfilter                        \
                libavcodec                         \
                libswresample                      \
                libswscale                         \
                libavutil                          \

CFLAGS += -Wall -g
CFLAGS := $(shell pkg-config --cflags $(FFMPEG_LIBS)) $(CFLAGS)
LDLIBS := $(shell pkg-config --libs $(FFMPEG_LIBS)) $(LDLIBS)

EXAMPLES=       decode_video_test                  \
                edge_detector_test                 \
                shape_preserving_test              \
                mv_estimation_test                 \

OBJS=$(addsuffix .o,$(EXAMPLES))

# the following examples make explicit use of the math library
edge_detector_test:           LDLIBS += -lm
shape_preserving_test:        LDLIBS += -lm
mv_estimation_test:           LDLIBS += -lm
.phony: all clean-test clean

all: $(OBJS) $(EXAMPLES)

clean-test:
	$(RM) test*.pgm test.h264 test.mp2 test.sw test.mpg

clean: clean-test
	$(RM) $(EXAMPLES) $(OBJS)
