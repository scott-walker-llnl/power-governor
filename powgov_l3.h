#pragma once
#include <stdint.h>
#include "powgov.h"
#include "powgov_profiles.h"

#define MAX_L3_SEQ 32

struct powgov_runtime;

struct powgov_l3
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
	uint64_t baseline_ipc;
	double scalability;
	short seq_end;
	struct phase_profile sequence[MAX_L3_SEQ];
	// unsigned char graph[MAX_PROFILES][MAX_PROFILES];
};

void l3_analysis(struct powgov_runtime *runtime);
