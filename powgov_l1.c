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
	this_profile.frq = ((double) (runtime->sampler->l1->new_sample.aperf -
				runtime->sampler->l1->prev_sample.aperf) /
				(double) (runtime->sampler->l1->new_sample.mperf -
				runtime->sampler->l1->prev_sample.mperf)) * runtime->sys->max_non_turbo;
	this_profile.occurrences = 1;

	// TODO: can we not do this every iteration?
	update_max(runtime, &this_profile);

	// if there is not a current workload, we should use this profile as the current workload
	// TODO: we need to initialize current workload to zero in main
	if (runtime->sampler->l1->current_workload.ipc == 0.0)
	{
		printf("INITIALIZE CURRENT\n");
		runtime->sampler->l1->current_workload = this_profile;
		runtime->sampler->l1->current_workload.occurrences = 1;
		runtime->sampler->l1->phase_begin = runtime->sampler->l1->new_sample;
		runtime->sampler->l1->current_workload.frq =
				(double) ((runtime->sampler->l1->new_sample.frq_data>> 8) & 0xFFul);
		print_profile(&runtime->sampler->l1->current_workload);
		// return here, this is first iteration so we don't have enough info yet
		return;
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
	classify_workload(runtime, &this_profile);

	// check control flow, did we stay in the same phase or has it changed?
	// TODO: check if we are in L3 predicted phase next
	if (branch_same_workload(runtime, &this_profile))
	{
		// classify current workload and react accordingly
		if (runtime->cfg->man_cpu_ctrl)
		{
			react_to_workload(runtime);
		}
	}
	else
	{
		branch_change_workload(runtime, &this_profile);
		if (runtime->cfg->man_cpu_ctrl)
		{
			react_to_workload(runtime);
		}
	}
}

void update_current_workload(struct powgov_runtime *runtime, struct workload_profile * prof)
{
	printf("\nupdating workload profile\nnew\n");
	print_profile(prof);
	struct workload_profile *current = &runtime->sampler->l1->current_workload;
	printf("pre\n");
	print_profile(current);
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
	printf("post\n");
	print_profile(current);
}

int branch_same_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile)
{
	struct workload_profile scaled_profile;
	double dist_to_recent;

	// scale the just seen workload to the recently executed workload
	frequency_scale_phase(this_profile, this_profile->frq,
			runtime->sampler->l1->current_workload.frq, &scaled_profile);
	// calculate distance between epoch and last seen phase
	dist_to_recent = workload_metric_distance(&scaled_profile,
			&runtime->sampler->l1->current_workload, &runtime->classifier->prof_maximums);

#ifdef NEW_DEBUG
	printf("\ndistance to recent %lf\ncurrent\n", dist_to_recent);
	print_profile(&runtime->sampler->l1->current_workload);
	printf("new\n");
	print_profile(this_profile);
#endif

	if (dist_to_recent < runtime->classifier->dist_thresh)
	{
		// we are in the same phase, update it
		update_current_workload(runtime, &scaled_profile);
		return 1;
	}
	// we are not in the same phase
	return 0;
}

int branch_change_workload(struct powgov_runtime *runtime, struct workload_profile *this_profile)
{
	// we are in a new phase, we can now say what the last phase was 
	struct phase_profile this_phase;
	this_phase.cycles = runtime->sampler->l1->prev_sample.tsc_data -
				runtime->sampler->l1->phase_begin.tsc_data;
	this_phase.workload = runtime->sampler->l1->current_workload;
	this_phase.phase_occurrences = 1;

	// update phase begin marker with the last sample, which was the end
	runtime->sampler->l1->phase_begin = runtime->sampler->l1->prev_sample;

	// check if current_workload has been predicted by L3 
	if (runtime->sampler->l3->predicted_phase != NULL)
	{
		struct workload_profile scaled_profile;
		double dist_to_predicted;
		// scale the epoch phase to the frequency of the last seen phase
		frequency_scale_phase(this_profile, this_profile->frq,
				runtime->sampler->l3->predicted_phase->workload.frq, &scaled_profile);
		// calculate distance between epoch and last seen phase
		dist_to_predicted = workload_metric_distance(&scaled_profile,
				&runtime->sampler->l3->predicted_phase->workload,
				&runtime->classifier->prof_maximums);
		if (dist_to_predicted < runtime->classifier->dist_thresh)
		{
			// prediction was correct
			runtime->sampler->l1->current_workload.frq_target =
					runtime->sampler->l3->predicted_phase->workload.frq_target;
		}
		else
		{
			// prediction was incorrect, we don't know what to do so start at defaults and
			// learn what is needed
			runtime->sampler->l1->current_workload.frq_target = runtime->sys->max_pstate;
		}
	}

	// a phase just completed and the current_workload represents that phase
	// now we need to find the best matching previously seen phase to update the registry
	struct phase_profile *phases = runtime->classifier->phases;
	struct workload_profile min_scaled_profile;
	double min_dist = DBL_MAX;
	int numunder = 0;
	int min_idx = -1;
	int i;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		// measure the distance to every known phase
		frequency_scale_phase(&runtime->sampler->l1->current_workload, 
				runtime->sampler->l1->current_workload.frq, phases[i].workload.frq,
				&this_phase.workload);
		double dist = phase_metric_distance(&this_phase, &phases[i],
				&runtime->classifier->prof_maximums, runtime->sampler->l3->minimum_cycles,
				runtime->sampler->l3->maximum_cycles);
#ifdef NEW_DEBUG
		printf("current phase dist from phase %d: %lf\n", i, dist);
		print_profile(&phases[i].workload);
		print_profile(&this_phase.workload);
#endif
		// if we are within the threshold then we have found a matching phase
		if (dist < min_dist)
		{
			min_dist = dist;
			min_idx = i;
			if (dist < runtime->classifier->dist_thresh)
			{
				numunder++;
			}
			min_scaled_profile = this_phase.workload;
		}
	}
	if (min_idx >= 0 && min_dist < runtime->classifier->dist_thresh)
	{
		// the phase that just completed HAS been seen before
		struct phase_profile *identified_phase = update_phase(runtime, &min_scaled_profile,
				&phases[min_idx], this_phase.cycles);
		update_graph_node(runtime, identified_phase);
		update_minmax_cycles(runtime, this_phase.cycles);
	}
	else
	{
		// the phase that just completed HAS NOT been seen before
		if (runtime->classifier->numphases >= MAX_PROFILES)
		{
			printf("ERROR: out of profile storage, increase the limit or change the sensitivity\n");
			return -1;
		}
		struct phase_profile *new_phase = add_phase(runtime, 
				&runtime->sampler->l1->current_workload, this_phase.cycles);
		printf("added new phase\n");
		print_profile(&new_phase->workload);
		add_graph_node(runtime, new_phase);
		update_minmax_cycles(runtime, this_phase.cycles);
	}
	printf("number of phases %d\n", runtime->classifier->numphases);

	// TODO: if there are many matches, we should combine similar phases
	/* if (numunder > 0) */
	/* { */
	/* 	agglomerate_profiles(runtime); */
	/* } */

	// the current workload has changed, replace it
	runtime->sampler->l1->current_workload = *this_profile;
	return 0;
}

void react_to_workload(struct powgov_runtime *runtime)
{
	struct workload_profile *profile = &runtime->sampler->l1->current_workload;

	// do frequency override if enabled
	if (runtime->cfg->mem_frq_override != 0.0 && profile->class == CLASS_MEM)
	{
		profile->frq_target = runtime->cfg->mem_frq_override;
		set_perf(runtime, (uint16_t) profile->frq_target);
		return;
	}
	else if (runtime->cfg->cpu_frq_override != 0.0 && profile->class == CLASS_CPU)
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
			profile->frq_target += runtime->cfg->frq_change_step;
			profile->unthrottle_cycles = 0;
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
	profile->unthrottle_cycles++;
}
