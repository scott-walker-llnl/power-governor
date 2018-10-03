#pragma once
#include "powgov_sampler.h"

struct powgov_l1
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
};


void pow_aware_perf(struct powgov_runtime *runtime);
