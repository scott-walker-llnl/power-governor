#pragma once
#include "powgov.h"
#define MAX_PROFILES 20
#define NUM_CLASSES 4

#define CLASS_CPU_SLOPE_IPC 0.63396
#define CLASS_CPU_SLOPE_EPC 0.13005
#define CLASS_MEM_SLOPE_IPC 0.07642
/* #define CLASS_MEM_SLOPE_EPC 0.73337 */
#define CLASS_MEM_SLOPE_EPC 0.79
#define CLASS_MIX_SLOPE_IPC ((CLASS_CPU_SLOPE_IPC + CLASS_MEM_SLOPE_IPC) / 2)
#define CLASS_MIX_SLOPE_EPC ((CLASS_CPU_SLOPE_EPC + CLASS_MEM_SLOPE_EPC) / 2 )
#define CLASS_CPU_INTERCEPT_IPC 0.10806
#define CLASS_CPU_INTERCEPT_EPC -0.12874
#define CLASS_MEM_INTERCEPT_IPC 0.08295
#define CLASS_MEM_INTERCEPT_EPC -0.20863
#define CLASS_MIX_INTERCEPT_IPC ((CLASS_CPU_INTERCEPT_IPC + CLASS_MEM_INTERCEPT_IPC) / 2)
#define CLASS_MIX_INTERCEPT_EPC ((CLASS_CPU_INTERCEPT_EPC + CLASS_MEM_INTERCEPT_EPC) / 2 )

enum CLASS_ID
{
	CLASS_CPU,
	CLASS_MEM,
	CLASS_IO,
	CLASS_MIX,
	CLASS_UNKNOWN,
};

struct powgov_runtime;

struct phase_profile
{
	// these are all AVERAGES
	double ipc; // instructions/cycle measured for this phase
	double mpc; // LLC miss/cycle measured for this phase
	double rpc; // resource stalls/cycle measured for this phase
	double epc; // execution stalls/cycle measured for this phase
	double bpc; // branch instructions/cycle measured for this phase

	double avg_frq; // average frequency measured for this phase
	uint16_t frq_high; // max frequency measured for this phase
	uint16_t frq_low; // min frequency measured for this phase
	double frq_target; // what frequency in *100MHz the algorithm thinks phase should run at
	double avg_cycle; // average number of cycles it takes to execute this phase
	uint32_t num_throttles; // number of times this phase was throttled last time (aka misprediction)
	uint64_t occurrences; // how many times this phase was detected
	char lastprev;
	char class;
	char unthrot_count;
	char reclass_count;
	char frq_duty_count;
};

struct powgov_classifier
{
	double dist_thresh; // the threshold for profile cluster identification
	double pct_thresh; // the threshold for displaying dumped profiles as execution percent
	int numphases; // the current number of phases
	struct phase_profile profiles[MAX_PROFILES]; // the cluster centers
	//float transition_table[MAX_PROFILES][MAX_PROFILES]; // phase transition matrix
	struct phase_profile prof_maximums; // minimum sampled values used for scaling and normalization
	struct phase_profile prof_minimums; // minimum sampled values used for scaling and normalization
	struct phase_profile prof_class[NUM_CLASSES]; // pre-computed values for various classes of workload
	int recentphase;
	double mem_freq_throttle;
};


double metric_distance(struct phase_profile *old, struct phase_profile *new, struct phase_profile *maximums, struct phase_profile *minimums);
void agglomerate_profiles(struct powgov_runtime *runtime);
void remove_unused(struct powgov_runtime *runtime);
void update_minmax(struct phase_profile *this_profile, struct phase_profile *maximums, 
		struct phase_profile *minimums);
void print_profile(struct phase_profile *prof);
void update_profile(struct powgov_runtime *runtime, struct phase_profile *this_profile, int profidx, uint64_t perf, unsigned this_throttle, double avgfrq, int lastphase);
void add_profile(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, unsigned this_throttle, double avgfrq, int lastphase);
int classify_phase(struct powgov_runtime *runtime, struct phase_profile *phase, uint64_t perf);
void classify_and_react(struct powgov_runtime *runtime, int phase, char wasthrottled, uint64_t perf);
double ipc_scale(double ipc_unscaled, double frq_source, double frq_target);
void frequency_scale_phase(struct phase_profile *unscaled_profile, double frq_source, double frq_target, struct phase_profile *scaled_profile);
int branch_same_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq);
int branch_change_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq);
