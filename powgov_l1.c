#include <float.h>
#include <stdint.h>
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"
#include "master.h"
#include "powgov_util.h"

void l1_analysis(struct powgov_runtime *runtime)
{
	// 1. find phase barriers
	// 2. avoid throttles (happens in "classify")
	// 3. set frequency (happens in "classify")
	// 4. check L3 prediction (happens in branch_change_workload)

	// TODO: this should actually be used...
	if (runtime->sampler->l1->new_sample.energy_data < 
			runtime->sampler->l1->prev_sample.energy_data)
	{
		runtime->power->energy_overflow++;
	}

	// gather the metrics for the last executed epoch
	uint64_t this_instret = runtime->sampler->l1->new_sample.instret -
		runtime->sampler->l1->prev_sample.instret;
	uint64_t this_cycle = runtime->sampler->l1->new_sample.tsc_data -
		runtime->sampler->l1->prev_sample.tsc_data;
	uint64_t this_llcmiss = runtime->sampler->l1->new_sample.llcmiss -
		runtime->sampler->l1->prev_sample.llcmiss;
	uint64_t this_restalls = runtime->sampler->l1->new_sample.restalls -
		runtime->sampler->l1->prev_sample.restalls;
	uint64_t this_exstalls = runtime->sampler->l1->new_sample.exstalls -
		runtime->sampler->l1->prev_sample.exstalls;
	uint64_t this_branchret = runtime->sampler->l1->new_sample.branchret -
		runtime->sampler->l1->prev_sample.branchret;

	// create a profile for the last executed epoch
	struct workload_profile this_profile;
	this_profile.ipc = ((double) this_instret) / ((double) this_cycle);
	this_profile.mpc = ((double) this_llcmiss) / ((double) this_cycle);
	this_profile.rpc = ((double) this_restalls) / ((double) this_cycle);
	this_profile.epc = ((double) this_exstalls) / ((double) this_cycle);
	this_profile.bpc = ((double) this_branchret) / ((double) this_cycle);
	this_profile.avgfrq = ((double) (runtime->sampler->l1->new_sample.aperf -
				runtime->sampler->l1->prev_sample.aperf) /
				(double) (runtime->sampler->l1->new_sample.mperf -
				runtime->sampler->l1->prev_sample.mperf)) * runtime->sys->max_non_turbo;

	// TODO: can we not do this every iteration?
	update_minmax(runtime, &this_profile);

	// if there is not a current workload, we should use this profile as the current workload
	// TODO: we need to initialize current workload to zero in main
	if (runtime->sampler->l1->current_workload.ipc == 0.0)
	{
		runtime->sampler->l1->current_workload = this_profile;
		runtime->sampler->l1->current_workload.occurrences++;
	}

	// check if the last epoch was throttled by RAPL
	uint64_t limreasons = runtime->sampler->l1->new_sample.perflimit;
	if (limreasons & LIMIT_LOG_RAPL)
	{
		// clear out the throttle flags, so they will indicate again next time
		write_msr_by_coord(0, 0, 0, MSR_CORE_PERF_LIMIT_REASONS, (limreasons & LIMIT_LOG_MASK));
		runtime->sampler->l1->isthrottled = 1;
	}

	// TODO avoid reclassifying every timestep
	int class = classify_workload(runtime, profile);

	// check control flow, did we stay in the same phase or has it changed?
	// TODO: check if we are in L3 predicted phase next
	if (branch_same_workload(runtime, &this_profile))
	{
		// classify current workload and react accordingly
		if (runtime->cfg->man_cpu_ctrl)
		{
			react_to_workload(runtime, &runtime->sampler->l1->current_workload);
		}
	}
	else
	{
		branch_change_workload(runtime, &this_profile);
		if (runtime->cfg->man_cpu_ctrl)
		{
			react_to_workload(runtime, &doop);
		}
	}
}

void update_current_workload(struct powgov_runtime *runtime, struct workload_profile * prof)
{
	struct workload_profile *current = &runtime->sampler->l1->current_workload;
	current->occurrences++;
	current->ipc = (current->ipc * (current->occurrences - 1.0 / current->occurrences) +
			prof->ipc * (1.0 / current->occurrences));
	current->mpc = (current->mpc * (current->occurrences - 1.0 / current->occurrences) +
			prof->mpc * (1.0 / current->occurrences));
	current->rpc = (current->rpc * (current->occurrences - 1.0 / current->occurrences) +
			prof->rpc * (1.0 / current->occurrences));
	current->epc = (current->epc * (current->occurrences - 1.0 / current->occurrences) +
			prof->epc * (1.0 / current->occurrences));
	current->bpc = (current->bpc * (current->occurrences - 1.0 / current->occurrences) +
			prof->bpc * (1.0 / current->occurrences));
	current->frq = (current->frq * (current->occurrences - 1.0 / current->occurrences) +
			prof->frq * (1.0 / current->occurrences));
}

int branch_same_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile)
{
	struct workload_profile scaled_profile;
	double dist_to_recent;

	// scale the epoch phase to the frequency of the last seen phase
	frequency_scale_phase(this_profile, this_profile->avgfrq,
			runtime->sampler->l1->current_workload.avg_frq, &scaled_profile);
	// calculate distance between epoch and last seen phase
	dist_to_recent = workload_metric_distance(&scaled_profile,
			runtime->sampler->l1->current_workload, &runtime->classifier->prof_maximums,
			&runtime->classifier->prof_minimums);

	if (dist_to_recent < runtime->classifier->dist_thresh)
	{
		// we are in the same phase, update it
		update_current_workload(runtime, &scaled_profile);
		// TODO: do we really need to re-classify if phase hasn't changed?
		return 1;
	}
	// we are not in the same phase
	return 0;
}

int branch_change_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile)
{
	// we are in a new phase, report previous phase info to L3
	double last_phase_cycles = runtime->sampler->l1->prev_sample.tsc_data -
				runtime->sampler->l1->phase_begin.tsc_data;
	if (runtime->sampler->l3->seq_end < MAX_L3_SEQ - 1)
	{
		runtime->sampler->l3->seq_end++;
		struct phase_profile *phase = 
				&runtime->sampler->l3->sequence[runtime->sampler->l3->seq_end];
		phase->workload = runtime->sampler->l1->current_workload;
		phase->cycles = last_phase_cycles;
		runtime->sampler->l1->phase_begin = runtime->sampler->l1->prev_sample;
	}

	// check if current_workload has been predicted by L3 
	struct workload_profile scaled_profile;
	double dist_to_predicted;
	// scale the epoch phase to the frequency of the last seen phase
	frequency_scale_phase(this_profile, this_profile->avgfrq,
			runtime->sampler->l3->predicted_phase.avg_frq, &scaled_profile);
	// calculate distance between epoch and last seen phase
	dist_to_predicted = workload_metric_distance(&scaled_profile,
			runtime->sampler->l3->predicted_phase.workload, &runtime->classifier->prof_maximums,
			&runtime->classifier->prof_minimums);
	// TODO: may not want to classify based on threshold for this
	if (dist_to_predicted < runtime->classifier->dist_thresh)
	{
		// prediction was correct
		runtime->sampler->current_workload.frq_target =
				runtime->sampler->l3->predicted_phase.workload.frq_target;
	}
	else
	{
		// prediction was incorrect, we don't know what to do so start at defaults and
		// learn what is needed
		runtime->sampler->current_workload.frq_target = runtime->sys->max_pstate;
	}

	// use last seen workload to update whichever phase it belongs to
	// search for the best match
	struct phase_profile *phases = runtime->classifier->phases;
	struct workload_profile min_scaled_profile;
	double min_dist = DBL_MAX;
	int numunder = 0;
	int min_idx = -1;
	int i;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		// measure the distance to every known phase
		struct workload_profile scaled_profile;
		frequency_scale_phase(runtime->sampler->l1->current_workload, 
				runtime->sampler->l1->current_workload.avgfrq, phases[i].workload.avg_frq,
				&scaled_profile);
		// TODO: this should use a separate metric distance that accounts for phase length
		double dist = phase_metric_distance(&scaled_profile, &phases[i].workload,
				&runtime->classifier->prof_maximums, &runtime->classifier->prof_minimums);
		// if we are within the threshold then we have found a matching phase
		if (dist < min_dist)
		{
			min_dist = dist;
			min_idx = i;
			if (dist < runtime->classifier->dist_thresh)
			{
				numunder++;
			}
			min_scaled_profile = scaled_profile;
		}
	}
	if (min_idx >= 0 && min_dist < runtime->classifier->dist_thresh)
	{
		// we found an existing phase which matches the phase which just completed 
		update_phase(runtime, &min_scaled_profile, &phases[min_idx], last_phase_cycles);
	}
	else
	{
		// the currently executing workload has never been seen before
		if (runtime->classifier->numphases >= MAX_PROFILES)
		{
			printf("ERROR: out of profile storage, increase the limit or change the sensitivity\n");
			return -1;
		}
		add_phase(runtime, runtime->sampler->l1->current_workload, last_phase_cycles);
	}

	// TODO: if there are many matches, we should combine similar phases
	/* if (numunder > 0) */
	/* { */
	/* 	agglomerate_profiles(runtime); */
	/* } */

	// the current workload has changed, replace it
	runtime->sampler->l1->current_workload = this_profile;
	return 0;
}

void react_to_workload(struct powgov_runtime *runtime, struct workload_profile *profile)
{
	
	struct phase_profile *phase = runtime->sampler->l3->current_phase;

	// do frequency override if enabled
	if (runtime->cfg->mem_frq_override != 0.0 && class == CLASS_MEM)
	{
		profile->frq_target = runtime->cfg->mem_frq_override;
		set_perf(runtime, (uint16_t) profile->frq_target);
		return;
	}
	else if (runtime->cfg->cpu_frq_override != 0.0 && class == CLASS_CPU)
	{
		profile->frq_target = runtime->cfg->cpu_frq_override;
		set_perf(runtime, (uint16_t) profile->frq_target);
		return;
	}

	// do throttle avoidance if enabled
	if (runtime->cfg->throttle_avoid && (runtime->sampler->l1->isthrottled ||
			runtime->sampler->l2->excursion))
	{
		// we were throttled use multiplicative decrease to drop frequency fast
		profile->frq_target -= runtime->cfg->frq_change_step *
				runtime->cfg->throttle_step_mult;
		profile->unthrottle_cycles = 0;
	}
	else
	{
		// we were not throttled, use additive increase to raise frequency slowly
		if (profile->frq_target < runtime->sys->max_pstate && profile->unthrottle_cycles >=
				runtime->cfg->fup_timeout)
		{
			profile->frq_target += FRQ_CHANGE_STEP;
			profiles[phase].unthrot_count = 0;
		}
	}

	// set the frequency
	// frequency duty cycling for fractional frequencies
	double ratio = profile->frq_target - (double)((uint16_t) profile->frq_target);
	if (profile->frq_duty_count < ratio * runtime->cfg->frq_duty_length)
	{
		// setting to floor first is safer for throttle avoidance
		// set perf to floor
		set_perf(runtime, (uint16_t) profile->frq_target);
	}
	else
	{
		// set perf to ciel
		set_perf(runtime, ((uint16_t) profile->frq_target) + 1);
	}
	profile->frq_duty_count = profile->frq_duty_count % runtime->cfg->frq_duty_length; 
	profile->frq_duty_count++;
	profile->unthrot_count++;
}
