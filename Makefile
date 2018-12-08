CC=g++

CPPFLAGS=-Wall -Wno-unused-label -O3 -g
LDFLAGS=-g -lm

SRCS=main_pc.cpp gauss_var.cpp gauss_fixed.cpp bench_common.cpp io.cpp bench_pc.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
HEADERS=bench_common.h gauss.h bench.h gauss_var.h gauss_fixed.h

all:gauss

# generic .o realisation
%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

gauss: $(OBJS)
	echo $(OBJS)
	$(CC) $(LDFLAGS) -o gauss $(OBJS)

clean:
	rm gauss -f
	rm $(OBJS) -f

make bench: clean gauss
	ls -lh gauss
	./gauss

# b:
# 	mkdir -p b