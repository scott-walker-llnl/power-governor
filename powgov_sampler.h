#pragma once
#include "powgov.h"

struct powgov_runtime;
struct powgov_l1;
struct powgov_l2;
struct powgov_l3;

struct data_sample
{
	uint64_t frq_data;
	uint64_t tsc_data;
	uint64_t aperf;
	uint64_t mperf;
	uint64_t energy_data;
	// uint64_t rapl_throttled;
	uint64_t therm;
	uint64_t perflimit;
	uint64_t instret;
	uint64_t llcmiss;
	uint64_t restalls;
	uint64_t exstalls;
	uint64_t branchret;
};

struct powgov_sampler
{
	unsigned long *samplectrs; // counter for number of samples on each thread
	struct data_sample **thread_samples;
	unsigned long numsamples;
	struct data_sample first_sample;
	struct powgov_l1 *l1;
	struct powgov_l2 *l2;
	struct powgov_l3 *l3;
	unsigned sps;
	unsigned long total_samples;
};

int sample_data(struct powgov_runtime *runtime);
void init_sampling(struct powgov_runtime *runtime);
