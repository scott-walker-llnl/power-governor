#include <float.h>
#include <stdint.h>
#include "powgov_l1.h"
#include "master.h"

void pow_aware_perf(struct powgov_runtime *runtime)
{
	static uint64_t begin_aperf = 0, begin_mperf = 0;
	static uint64_t last_aperf = 0, last_mperf = 0;
	static uint64_t tsc_timer = 0;
	/* static uint64_t phase_start_tsc = 0; */
	/* struct phase_profile *profiles = runtime->classifier->profiles; */

	if (last_aperf == 0)
	{
		read_msr_by_coord(0, 0, 0, MSR_IA32_MPERF, &begin_mperf);
		read_msr_by_coord(0, 0, 0, MSR_IA32_APERF, &begin_aperf);
		read_msr_by_coord(0, 0, 0, IA32_TIME_STAMP_COUNTER, &tsc_timer);
		/* phase_start_tsc = tsc_timer; */
		last_aperf = begin_aperf;
		last_mperf = begin_mperf;
	}

	uint64_t aperf, mperf;
	read_msr_by_coord(0, 0, 0, MSR_IA32_MPERF, &mperf);
	read_msr_by_coord(0, 0, 0, MSR_IA32_APERF, &aperf);

	if (runtime->sampler->l1->new_sample.energy_data < runtime->sampler->l1->prev_sample.energy_data)
	{
		runtime->power->energy_overflow++;
	}

	uint64_t this_instret = runtime->sampler->l1->new_sample.instret - runtime->sampler->l1->prev_sample.instret;
	uint64_t this_cycle = runtime->sampler->l1->new_sample.tsc_data - runtime->sampler->l1->prev_sample.tsc_data;
	//unsigned this_throttle = (runtime->sampler->l1->new_sample.rapl_throttled & 0xFFFFFFFF) -
	//		(runtime->sampler->l1->prev_sample.rapl_throttled & 0xFFFFFFFF);
	uint64_t this_llcmiss = runtime->sampler->l1->new_sample.llcmiss - runtime->sampler->l1->prev_sample.llcmiss;
	uint64_t this_restalls = runtime->sampler->l1->new_sample.restalls - runtime->sampler->l1->prev_sample.restalls;
	uint64_t this_exstalls = runtime->sampler->l1->new_sample.exstalls - runtime->sampler->l1->prev_sample.exstalls;
	uint64_t this_branchret = runtime->sampler->l1->new_sample.branchret - runtime->sampler->l1->prev_sample.branchret;

	//double total_avgfrq = ((double) (aperf - begin_aperf) / (double) (mperf - begin_mperf)) * runtime->sys->max_non_turbo;
	double phase_avgfrq = ((double) (aperf - last_aperf) / (double) (mperf - last_mperf)) * runtime->sys->max_non_turbo;
	uint64_t perf = ((runtime->sampler->l1->new_sample.frq_data & 0xFFFFul) >> 8);

	struct phase_profile this_profile;

	this_profile.ipc = ((double) this_instret) / ((double) this_cycle); // instructions per cycle
	this_profile.mpc = ((double) this_llcmiss) / ((double) this_cycle); // cache misses per cycle
	this_profile.rpc = ((double) this_restalls) / ((double) this_cycle); // resource stalls per cycle
	this_profile.epc = ((double) this_exstalls) / ((double) this_cycle); // execution stalls per cycle
	this_profile.bpc = ((double) this_branchret) / ((double) this_cycle); // branch instructions retired per cycle

	update_minmax(&this_profile, &runtime->classifier->prof_maximums, 
			&runtime->classifier->prof_minimums);

	if (runtime->classifier->numphases >= MAX_PROFILES)
	{
		remove_unused(runtime);
	}

	if (runtime->classifier->numphases == 0)
	{
		add_profile(runtime, &this_profile, perf, 0, phase_avgfrq, runtime->classifier->recentphase);
		return;
	}
	// we do phase analysis
	// if current execution is similar to previously seen phase, update that phase
	// check to see if we are in the same phase
	if (runtime->classifier->recentphase > runtime->classifier->numphases)
	{
		//printf("recent phase no longer exists\n");
		runtime->classifier->recentphase = -1;
	}

	uint64_t limreasons = runtime->sampler->l1->new_sample.perflimit;
	char isthrottled = 0;
	char wasthrottled = 0;
	if (limreasons & LIMIT_LOG_RAPL)
	{
		write_msr_by_coord(0, 0, 0, MSR_CORE_PERF_LIMIT_REASONS, (limreasons & LIMIT_LOG_MASK));
		wasthrottled = 1;
	}
	if (limreasons & LIMIT_ON_RAPL)
	{
		isthrottled = 1;
	}

	// we may be in the same phase
	if (branch_same_phase(runtime, &this_profile, perf, wasthrottled, isthrottled, phase_avgfrq))
	{
		return;
	}
	// the phase changed
	/* profiles[runtime->classifier->recentphase].avg_cycle = (double) */
		/* (runtime->sampler->l1->new_sample.tsc_data - phase_start_tsc); */

	/* double pcyc = (double) (runtime->sampler->l1->new_sample.tsc_data - phase_start_tsc); */
	/* profiles[runtime->classifier->recentphase].avg_cycle = */
	/* 	weighted_avg_flt(profiles[runtime->classifier->recentphase].occurrences, */
	/* 	profiles[runtime->classifier->recentphase].avg_cycle, */
	/* 	pcyc); */
	/* printf("cycle %d len: %lf\n", runtime->classifier->recentphase,  */
			/* profiles[runtime->classifier->recentphase].avg_cycle); */

	/* phase_start_tsc = runtime->sampler->l1->new_sample.tsc_data; */
	last_aperf = aperf;
	last_mperf = mperf;
	branch_change_phase(runtime, &this_profile, perf, wasthrottled, isthrottled, phase_avgfrq);
}

int branch_same_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq)
{
	if (runtime->classifier->recentphase < 0)
	{
		return 0;
	}

	struct phase_profile scaled_profile = *this_profile;
	double dist_to_recent;
	struct phase_profile *profiles = runtime->classifier->profiles;
	frequency_scale_phase(this_profile, phase_avgfrq, profiles[runtime->classifier->recentphase].avg_frq, &scaled_profile);
	dist_to_recent = metric_distance(&scaled_profile, &profiles[runtime->classifier->recentphase], &runtime->classifier->prof_maximums, &runtime->classifier->prof_minimums);

	if (dist_to_recent < runtime->classifier->dist_thresh)
	{
#ifdef DEBUG
		printf("Phase has not changed, dist to recent %lf\n", dist_to_recent);
#endif
		// we are in the same phase, update it
		update_profile(runtime, &scaled_profile, runtime->classifier->recentphase, perf, wasthrottled, phase_avgfrq, runtime->classifier->recentphase);
		classify_and_react(runtime, runtime->classifier->recentphase, wasthrottled, perf);
		return 1;
	}
#ifdef DEBUG
	printf("Phase has changed, dist to recent %lf\n", dist_to_recent);
	printf("\trecent phase id %d: ipc %lf, epc %lf, freq %lf. new phase: ipc %lf, epc %lf, freq %lf\n",
		runtime->classifier->recentphase, profiles[runtime->classifier->recentphase].ipc, profiles[runtime->classifier->recentphase].epc, profiles[runtime->classifier->recentphase].avg_frq,
		this_profile->ipc, this_profile->epc, phase_avgfrq);
#endif
	return 0;
}

int branch_change_phase(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, char wasthrottled, char isthrottled, double phase_avgfrq)
{
	// we are in a different phase so search for the best match
	double distances[MAX_PROFILES];
	int k;
	struct phase_profile *profiles = runtime->classifier->profiles;
	for (k = 0; k < runtime->classifier->numphases; k++)
	{
		distances[k] = DBL_MAX;
	}

	if (this_profile->ipc == 0.0 || this_profile->epc == 0.0)
	{
		printf("ERROR: bad profile data\n");
	}

#ifdef DEBUG
		printf("Search for most similar phase\n");
#endif

	struct phase_profile min_scaled_profile;
	double min_dist = DBL_MAX;
	int numunder = 0;
	int min_idx = -1;
	int i;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		// measure the distance to every known phase
		if (i == runtime->classifier->recentphase)
		{
			continue;
		}
		struct phase_profile scaled_profile;
		frequency_scale_phase(this_profile, phase_avgfrq, profiles[i].avg_frq, &scaled_profile);
		distances[i] = metric_distance(&scaled_profile, &profiles[i], &runtime->classifier->prof_maximums, &runtime->classifier->prof_minimums);

		if (distances[i] < min_dist)
		{
			min_dist = distances[i];
			min_idx = i;
			min_scaled_profile = scaled_profile;
		}
		if (distances[i] < runtime->classifier->dist_thresh)
		{
			numunder++;
		}
	}

#ifdef DEBUG
	int q;
	for (q = 0; q < runtime->classifier->numphases; q++)
	{
		printf("distance from %d to %d is %lf\n", runtime->classifier->numphases, q, distances[q]);
	}
#endif

	if (min_idx >= 0 && min_dist < runtime->classifier->dist_thresh)
	{
		// we found an existing phase which matches the currently executing workload
#ifdef DEBUG
		printf("Found existing phase, with distance %lf\n", min_dist);
#endif
		update_profile(runtime, &min_scaled_profile, min_idx, perf, isthrottled, phase_avgfrq, runtime->classifier->recentphase);
		runtime->classifier->recentphase = min_idx;
		classify_and_react(runtime, runtime->classifier->recentphase, wasthrottled, perf);
	}
	else
	{
		// the currently executing workload has never been seen before
		if (runtime->classifier->numphases >= MAX_PROFILES)
		{
			printf("ERROR: out of profile storage, increase the limit or change the sensitivity\n");
			return -1;
		}
		// add new phase
		add_profile(runtime, this_profile, perf, isthrottled, phase_avgfrq, runtime->classifier->recentphase);
		runtime->classifier->recentphase = runtime->classifier->numphases - 1;
	}

	// if there are many matches, we should combine similar phases
	if (numunder > 0)
	{
		agglomerate_profiles(runtime);
	}
	return 0;
}
