##################################################
#
# Makefile for TV-Web
#
###################################################

CC = cc
FFTW_INC=-I/usr/include
FFTW_LIB=-L/usr/lib
CC_OPTIONS = -Wall $(FFTW_INC) $(FFTW_LIB)  -m64  #-DNEWFORMAT  #-DGADGET_UNITS #-DJUANK #-DNEWFORMAT #-DSWAP -m64
LNK_OPTIONS = -Wall -lm  -lfftw3 $(FFTW_LIB) $(FFTW_INC)
EXEC = smooth_bolshoi.x
EXEC_GLUE = glue.x
EXEC_ENV = get_env.x
EXEC_MIX = mix_env.x
LIBS = -lm

OBJS = smooth.o io.o endian.o field.o fft_diff.o tools.o READ_ART.o
INCL = Makefile.altix io.h field.h struct.h endian.h

OBJS_GLUE = glue_meshes.o io.o endian.o field.o fft_diff.o tools.o READ_ART.o
OBJS_ENV = get_env.o io.o endian.o field.o fft_diff.o tools.o READ_ART.o
OBJS_MIX = mix_env.o io.o endian.o field.o fft_diff.o tools.o READ_ART.o


.SUFFIXES: .o .or .c

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

all: $(EXEC)

$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) $(LNK_OPTIONS) -o $(EXEC)

glue: $(EXEC_GLUE)

$(EXEC_GLUE):$(OBJS_GLUE)
	$(CC) $(OBJS_GLUE) $(LIBS) $(LNK_OPTIONS) -o $(EXEC_GLUE)

env: $(EXEC_ENV)

$(EXEC_ENV):$(OBJS_ENV)
	$(CC) $(OBJS_ENV) $(LIBS) $(LNK_OPTIONS) -o $(EXEC_ENV)

mix: $(EXEC_MIX)

$(EXEC_MIX):$(OBJS_MIX)
	$(CC) $(OBJS_MIX) $(LIBS) $(LNK_OPTIONS) -o $(EXEC_MIX)




clean:
	rm -f $(OBJS) *~ core* ${EXEC}

