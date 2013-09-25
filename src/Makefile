PROG = debi 

OBJS = debi.o ItemsetOutput.o Bitmap.o BaseBitmap.o 

CC = g++ -g
CFLAGS = -ansi -W -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -O3  
LFLAGS = -lm 

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $^

%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@


clean:
	-rm -f $(OBJS) $(PROG)
