all: src/blockbuster.c
	gcc -pedantic src/blockbuster.c -o bin/blockbuster.x -lm -std=c99 -D_GNU_SOURCE

blockbuster: src/blockbuster.c
	gcc -pedantic src/blockbuster.c -o bin/blockbuster.x -lm -std=c99 -D_GNU_SOURCE

clean:
	rm blockbuster.x
