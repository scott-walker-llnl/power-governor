#pragma once

struct data_sample
{
	uint64_t frq_data;
	uint64_t tsc_data;
	uint64_t energy_data;
	uint64_t rapl_throttled;
	uint64_t therm;
	uint64_t perflimit;
	uint64_t instret;
	uint64_t llcmiss;
	uint64_t restalls;
	uint64_t exstalls;
	uint64_t branchret;
};

int sample_data(struct powgov_runtime *runtime);
void init_sampling(struct powgov_runtime *runtime);
