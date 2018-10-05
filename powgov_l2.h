#pragma once
#include "powgov.h"

struct powgov_runtime;

struct powgov_l2
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
	char excursion;
};

void l2_analysis(struct powgov_runtime *runtime);
