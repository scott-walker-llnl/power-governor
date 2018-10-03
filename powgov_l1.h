#pragma once
#include "powgov.h"

#define LIMIT_LOG_RAPL 0xC000000
#define LIMIT_LOG_MASK 0xF3FFFFFF
#define LIMIT_ON_RAPL 0x0000C00
#define MSR_CORE_PERF_LIMIT_REASONS 0x64F

struct powgov_runtime;

struct powgov_l1
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
};

void pow_aware_perf(struct powgov_runtime *runtime);
int branch_same_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq);
int branch_change_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq);
