#pragma once
#include "powgov_sampler.h"

struct powgov_l2
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
};

void l2_analysis(struct powgov_runtime *runtime);
