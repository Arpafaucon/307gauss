gauss: gauss.c bench.c io.c gauss.h main.h
	gcc gauss.c bench.c io.c -o gauss -Wall -O3 -g

clean:
	rm gauss