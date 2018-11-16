#pragma once
#include "powgov.h"

struct powgov_runtime;

struct powgov_l2
{
	struct data_sample new_sample; // the newly gathered sample at l2 sample rate
	struct data_sample prev_sample; // the last gathered sample at l2 rate
	unsigned long interval; // the sample rate of this level
	char excursion; // if true the power limit has been exceeded
};

void l2_analysis(struct powgov_runtime *runtime);
