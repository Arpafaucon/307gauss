CFLAGS = -Wall -Wno-unused-label -O3 -g
gauss: gauss.c bench.c io.c gauss.h main.h
	gcc gauss.c bench.c io.c -o gauss ${CFLAGS}

clean:
	rm gauss