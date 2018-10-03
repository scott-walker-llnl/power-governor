#include <string.h>
#include "powgov_l3.h"

void l3_analysis(struct powgov_runtime *runtime)
{
	// macro performance view
	// if counter > N seconds
	// then
	//     if exec trace is the same
	//     then
	//         if macro performance >= baseline
	//         then
	//             throttle is OK, keep throttling or throttle more
	//         else
	//             throttle is bad, unthrottle
	//     else
	//         we are not executing the same thing, new macro comparison
	// TODO if exec trace is the same (higher level classification?)
	// TODO only need l1->new_sample
	runtime->sampler->l3->prev_sample = runtime->sampler->l3->new_sample;
	runtime->sampler->l3->new_sample = runtime->sampler->l1->new_sample;
	
	printf("\nphase sequence: ");
	// use phase sequence to find if
	//     cpu and non-cpu phases alternate
	//     how often do they alternate
	//     want to find non-cpu phase before cpu phase
	//     how often does mem phase not lead into cpu phase
	int i;
	for (i = 0; i < runtime->sampler->l3->seq_end; i++)
	{
		printf("%d (%lu), ", runtime->sampler->l3->sequence[i], 
				runtime->sampler->l3->seq_cycles[i]);
		if (i > 0)
		{
			// TODO: agglomerate and remove needs to move graph columns
			/* runtime->sampler->l3->graph[runtime->sampler->l3->sequence[i-1]] */
					/* [runtime->sampler->l3->sequence[i]] = 1; */
		}
	}
	/* for (i = 0; i < MAX_PROFILES; i++) */
	/* { */
	/* 	printf("%d, ", runtime->sampler->l3->sequence[i]); */
	/* } */
	/* printf("\n"); */
	/* for (i = 0; i < MAX_PROFILES; i++) */
	/* { */
	/* 	int j; */
	/* 	printf("%d {, ", runtime->sampler->l3->sequence[i]); */
	/* 	for (j = 0; j < MAX_PROFILES; j++) */
	/* 	{ */
	/* 		printf("%d, ", runtime->sampler->l3->graph[runtime->sampler->l3->sequence[i]] */
	/* 				[runtime->sampler->l3->sequence[j]]); */
	/* 	} */
	/* 	printf("}\n"); */
	/* } */
	printf("\n");
	// set the first phase to the last seen, so next iteration knows what came before
	runtime->sampler->l3->sequence[0] = 
			runtime->sampler->l3->sequence[runtime->sampler->l3->seq_end];
	runtime->sampler->l3->seq_cycles[0] = 
			runtime->sampler->l3->seq_cycles[runtime->sampler->l3->seq_end];
	memset(runtime->sampler->l3->sequence + 1, -1, (MAX_L3_SEQ - 1) * sizeof(unsigned char));
	memset(runtime->sampler->l3->seq_cycles + 1, 0, (MAX_L3_SEQ - 1) * sizeof(uint64_t));
	runtime->sampler->l3->seq_end = 0;
}
