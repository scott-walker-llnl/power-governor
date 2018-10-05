#pragma once
#include "powgov.h"
#define MAX_PROFILES 20
#define NUM_CLASSES 4
#define SCALE_OUTLIER_THRESH_LOW 0.8
#define SCALE_OUTLIER_THRESH_HIGH 1.2
#define SCALE_THRESH 1.0 // 100MHz

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

struct powgov_runtime;

enum CLASS_ID
{
	CLASS_CPU,
	CLASS_MEM,
	CLASS_IO,
	CLASS_MIX,
	CLASS_UNKNOWN,
};

/*
struct workload_profile
{
	// these are all AVERAGES
	double ipc; // instructions/cycle measured for this workload
	double mpc; // LLC miss/cycle measured for this workload
	double rpc; // resource stalls/cycle measured for this workload
	double epc; // execution stalls/cycle measured for this workload
	double bpc; // branch instructions/cycle measured for this workload

	double avg_frq; // average frequency measured for this workload
	uint16_t frq_high; // max frequency measured for this workload
	uint16_t frq_low; // min frequency measured for this workload
	double frq_target; // what frequency in *100MHz the algorithm thinks workload should run at
	double avg_cycle; // average number of cycles it takes to execute this workload
	uint32_t num_throttles; // number of times this workload was throttled last time (aka misprediction)
	uint64_t occurrences; // how many times this workload was detected
	uint64_t phase_occurrences; // how many times this phase occurred
	char lastprev;
	char class;
	char unthrot_count;
	char reclass_count;
	char frq_duty_count;
};
*/

struct workload_profile
{
	// these are all AVERAGES
	double ipc; // instructions/cycle measured for this workload
	double mpc; // LLC miss/cycle measured for this workload
	double rpc; // resource stalls/cycle measured for this workload
	double epc; // execution stalls/cycle measured for this workload
	double bpc; // branch instructions/cycle measured for this workload
	double frq; // average frequency measured for this workload
	double frq_target; // the frequency being used for this workload
	uint64_t occurrences; // how many times this workload has occurred sequentially
	char class; // what class is the current workload (cpu, mem, etc)
	char unthrottle_cycles; // how many cycles has the workload gone unthrottled
	char frq_duty_count; // used to duty cycle the processor frequency to 10's of MHz
};

struct phase_profile
{
	struct workload_profile workload;
	double cycles;
	unsigned short phase_occurrences;
};

struct powgov_classifier
{
	double dist_thresh; // the threshold for profile cluster identification
	double pct_thresh; // the threshold for displaying dumped profiles as execution percent
	int numphases; // the current number of phases
	struct phase_profile phases[MAX_PROFILES]; // the cluster centers
	//float transition_table[MAX_PROFILES][MAX_PROFILES]; // phase transition matrix
	struct workload_profile prof_maximums; // minimum sampled values used for scaling and normalization
	struct workload_profile prof_minimums; // minimum sampled values used for scaling and normalization
	struct workload_profile prof_class[NUM_CLASSES]; // pre-computed values for various classes of workload
	int recentphase;
	double mem_freq_throttle;
};


double workload_metric_distance(struct workload_profile *old, struct workload_profile *new, struct workload_profile *maximums, struct workload_profile *minimums);
void agglomerate_profiles(struct powgov_runtime *runtime);
void remove_unused(struct powgov_runtime *runtime);
void update_minmax(struct workload_profile *this_profile, struct workload_profile *maximums, 
		struct workload_profile *minimums);
void print_profile(struct workload_profile *prof);
void update_profile(struct powgov_runtime *runtime, struct workload_profile *this_profile, struct workload_profile *prof, uint64_t perf, double avgfrq);
void add_profile(struct powgov_runtime *runtime, struct workload_profile *this_profile, uint64_t perf, double avgfrq);
int classify_workload(struct powgov_runtime *runtime, struct workload_profile *phase, uint64_t perf);
void frequency_scale_phase(struct workload_profile *unscaled_profile, double frq_source, double frq_target, struct workload_profile *scaled_profile);
