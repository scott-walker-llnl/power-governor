INCLUDES_PATH=/home/walker8/local/include
LIBS_PATH=/home/walker8/local/lib
CFLAGS=-O3

powgov: powgov.c powgov.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

l1: powgov_l1.c powgov_l1.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

l2: powgov_l2.c powgov_l2.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

l3: powgov_l3.c powgov_l3.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

profiles: powgov_profiles.c powgov_profiles.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

sampler: powgov_sampler.c powgov_sampler.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

util: powgov_util.c powgov_util.h
	gcc -c $< ${CFLAGS} -o $@.o -I${INCLUDES_PATH} -std=c99 -Wall

all: powgov l1 l2 l3 profiles sampler util
	gcc $(addsuffix .o, $^) ${CFLAGS} -lm -L${LIBS_PATH} -lmsr -I${INCLUDES_PATH} -o power_governor -Wall

clean:
	/bin/rm -rf *.o
	/bin/rm -rf power_governor 
