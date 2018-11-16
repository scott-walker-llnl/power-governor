#pragma once
#include "powgov.h"
#include "powgov_profiles.h"

#define LIMIT_LOG_RAPL 0xC000000
#define LIMIT_LOG_MASK 0xF3FFFFFF
#define LIMIT_ON_RAPL 0x0000C00
#define RECLASSIFY_INTERVAL 20
#define MSR_CORE_PERF_LIMIT_REASONS 0x64F

struct powgov_runtime;
struct workload_profile;

struct powgov_l1
{
	char isthrottled; // true if the throttle reasons register is indicating RAPL
	struct data_sample new_sample; // the newest gathered sample
	struct data_sample prev_sample; // the sample from last poll
	struct data_sample phase_begin; // the sample from the beginning of this phase
	unsigned long interval; // the sample rate of this level
	struct workload_profile current_workload; // the current workload descriptor
};

void l1_analysis(struct powgov_runtime *runtime);
void update_current_workload(struct powgov_runtime *runtime, struct workload_profile * prof);
int branch_same_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile);
int branch_change_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile);
void react_to_workload(struct powgov_runtime *runtime);
struct phase_profile *update_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, struct phase_profile *prof, double phase_cycles);
struct phase_profile *add_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, double phase_cycles);
