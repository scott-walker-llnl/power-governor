INCLUDES_PATH=/home/walker8/local/include
LIBS_PATH=/home/walker8/local/lib

powgov: powgov.c powgov.h
	gcc -c $< -O3 -o $@.o -std=c99 -Wall

util: powgov_util.c powgov_util.h
	gcc -c $M -O3 -o $@.o -std=c99 -Wall

sampler: sampler.c
	gcc sampler.c -O3 -lm -L${LIBS_PATH} -lmsr -I${INCLUDES_PATH} -o sampler -std=c99 -Wunused

clean:
	/bin/rm -rf powgov.o
	/bin/rm -rf powgov.o
