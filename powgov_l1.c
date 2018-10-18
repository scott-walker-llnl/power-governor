#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
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

	struct data_sample *new_sample = &runtime->sampler->l1->new_sample;
	struct data_sample *prev_sample = &runtime->sampler->l1->prev_sample;

	uint64_t this_instret = new_sample->instret - prev_sample->instret;
	uint64_t this_cycle = new_sample->tsc_data - prev_sample->tsc_data;
	uint64_t this_llcmiss = new_sample->llcmiss - prev_sample->llcmiss;
	uint64_t this_restalls = new_sample->restalls - prev_sample->restalls;
	uint64_t this_exstalls = new_sample->exstalls - prev_sample->exstalls;
	uint64_t this_branchret = new_sample->branchret - prev_sample->branchret;

		// create a profile for the last executed epoch
	struct workload_profile this_profile;
	this_profile.ipc = ((double) this_instret) / ((double) this_cycle);
	this_profile.mpc = ((double) this_llcmiss) / ((double) this_cycle);
	this_profile.rpc = ((double) this_restalls) / ((double) this_cycle);
	this_profile.epc = ((double) this_exstalls) / ((double) this_cycle);
	this_profile.bpc = ((double) this_branchret) / ((double) this_cycle);
	this_profile.frq = ((double) (new_sample->aperf - prev_sample->aperf) /
						(double) (new_sample->mperf - prev_sample->mperf)) *
						runtime->sys->max_non_turbo;
	this_profile.occurrences = 1;
	this_profile.class = CLASS_UNKNOWN;

	/* printf("perf status was %lx\n", runtime->sampler->l1->new_sample.frq_data); */

	// TODO: can we not do this every iteration?
	update_max(runtime, &this_profile);
	// TODO avoid reclassifying every timestep
	classify_workload(runtime, &this_profile);

	// if there is not a current workload, we should use this profile as the current workload
	// TODO: we need to initialize current workload to zero in main
	if (runtime->sampler->l1->current_workload.ipc == 0.0)
	{
		printf("INITIALIZE CURRENT\n");
		runtime->sampler->l1->current_workload = this_profile;
		runtime->sampler->l1->current_workload.occurrences = 1;
		runtime->sampler->l1->phase_begin = runtime->sampler->l1->new_sample;
		runtime->sampler->l1->current_workload.frq =
				(double) ((runtime->sampler->l1->new_sample.frq_data >> 8) & 0xFFul);
		runtime->sampler->l1->current_workload.frq_target = runtime->sys->max_pstate;
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
	struct workload_profile *current = &runtime->sampler->l1->current_workload;
	current->occurrences++;

	struct data_sample *new_sample = &runtime->sampler->l1->new_sample;
	struct data_sample *phase_begin = &runtime->sampler->l1->phase_begin;
	uint64_t this_instret = new_sample->instret - phase_begin->instret;
	uint64_t this_cycle = new_sample->tsc_data - phase_begin->tsc_data;
	uint64_t this_llcmiss = new_sample->llcmiss - phase_begin->llcmiss;
	uint64_t this_restalls = new_sample->restalls - phase_begin->restalls;
	uint64_t this_exstalls = new_sample->exstalls - phase_begin->exstalls;
	uint64_t this_branchret = new_sample->branchret - phase_begin->branchret;

	// create a profile for the last executed epoch
	current->ipc = ((double) this_instret) / ((double) this_cycle);
	current->mpc = ((double) this_llcmiss) / ((double) this_cycle);
	current->rpc = ((double) this_restalls) / ((double) this_cycle);
	current->epc = ((double) this_exstalls) / ((double) this_cycle);
	current->bpc = ((double) this_branchret) / ((double) this_cycle);
	current->frq = ((double) (new_sample->aperf - phase_begin->aperf) /
					(double) (new_sample->mperf - phase_begin->mperf)) *
					runtime->sys->max_non_turbo;
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

	/* printf("\ndistance to recent %lf\n", dist_to_recent); */
#ifdef NEW_DEBUG
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
	// TODO: this_phase object confusing, it's workload data keeps being changed and
	//     its only use later is this_phase.cycles used to update/add new phase
	struct phase_profile this_phase;
	this_phase.cycles = runtime->sampler->l1->prev_sample.tsc_data -
				runtime->sampler->l1->phase_begin.tsc_data;
	//this_phase.workload = runtime->sampler->l1->current_workload;// write without read
	this_phase.phase_occurrences = 1;

	// update phase begin marker with the last sample, which was the end
	runtime->sampler->l1->phase_begin = runtime->sampler->l1->prev_sample;

	// we don't know what this phase is yet, take a conservative approach first
	this_profile->frq_target = runtime->sys->max_pstate;

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
			this_profile->frq_target =
					runtime->sampler->l3->predicted_phase->workload.frq_target;
		}
	}

	/* printf("new phase cycles %lf\t", this_phase.cycles); */

#ifdef NEW_DEBUG
	printf("\nCompare new phase: (%lf cyc)\n", this_phase.cycles);
	print_profile(&runtime->sampler->l1->current_workload);
#endif
	// a phase just completed and the current_workload represents that phase
	// now we need to find the best matching previously seen phase to update the registry
	// find nearest neighbor WORKLOAD, anything within cycle thresh is valid
	struct phase_profile *phases = runtime->classifier->phases;
	struct workload_profile min_scaled_profile;
	double min_dist = DBL_MAX;
	double min_cycles = DBL_MAX;
	int min_idx = -1;
	int i;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		double cycdiff = fabs(phases[i].cycles - this_phase.cycles);
		if (cycdiff < runtime->classifier->phase_cycle_thresh)
		{
			// measure the distance to every known phase
			frequency_scale_phase(&runtime->sampler->l1->current_workload, 
					runtime->sampler->l1->current_workload.frq, phases[i].workload.frq,
					&this_phase.workload);
			double dist = workload_metric_distance(&this_phase.workload, &phases[i].workload,
					&runtime->classifier->prof_maximums);

#ifdef NEW_DEBUG
			printf("current workload dist from phase %d: %lf (%lf cyc)\n", i, dist, cycdiff);
			print_profile(&phases[i].workload);
			print_profile(&this_phase.workload);
#endif
			// if we are within the threshold then we have found a matching phase
			if (dist < min_dist)
			{
				// we found workload that is similar, see if this phase cyc len is also similar
				min_dist = dist;
				min_cycles = cycdiff;
				min_idx = i;
				min_scaled_profile = this_phase.workload;
			}
		}
	}
	/* printf("min diff %lf\t", min_dist); */
#ifdef NEW_DEBUG
	printf("found min %d with dist %lf\n", min_idx, min_dist);
#endif
	if (min_idx >= 0 && min_dist <= runtime->classifier->dist_thresh)
	{
		/* printf("is seen\n"); */
#ifdef NEW_DEBUG
		print_profile(&phases[min_idx].workload);
#endif
		// the phase that just completed HAS been seen before
		struct phase_profile *identified_phase = update_phase(runtime, &min_scaled_profile,
				&phases[min_idx], this_phase.cycles);
		update_graph_node(runtime, identified_phase);
	}
	else
	{
		/*
		printf("\nadding new phase\n");
		print_profile(&this_phase.workload);
		printf("\tcycles%lf\n", this_phase.cycles);
		printf("most similar phase\n");
		print_profile(&phases[min_idx].workload);
		printf("cycle diff %lf, distance %lf\n", min_cycles, min_dist);
		for (i = 0; i < runtime->classifier->numphases; i++)
		{
			print_profile(&phases[i].workload);
			printf("\tcycles%lf\n", phases[i].cycles);
			printf("\tcycle diff%lf\n", fabs(phases[i].cycles - this_phase.cycles));
		}
		*/
#ifdef NEW_DEBUG
		printf("\nList Distances:\n");
		for (i = 0; i < runtime->classifier->numphases; i++)
		{
			double cycdiff = fabs(phases[i].cycles - this_phase.cycles);
			frequency_scale_phase(&runtime->sampler->l1->current_workload, 
					runtime->sampler->l1->current_workload.frq, phases[i].workload.frq,
					&this_phase.workload);
			double dist = workload_metric_distance(&this_phase.workload, &phases[i].workload,
					&runtime->classifier->prof_maximums);
			printf("\tcycdiff %lf, dist %lf\n", cycdiff, dist);
		}
#endif

		// the phase that just completed HAS NOT been seen before
		if (runtime->classifier->numphases >= MAX_PROFILES)
		{
			// try to eliminate unrepeated workloads
			//remove_unused(runtime);
			// there were no unrepeated workloads
			if (runtime->classifier->numphases >= MAX_PROFILES)
			{
				printf("ERROR: out of profile storage, increase the limit or change the sensitivity\n");
				dump_error_report(runtime);
				exit(-1);
				return -1;
			}
		}
		struct phase_profile *new_phase = add_phase(runtime, 
				&runtime->sampler->l1->current_workload, this_phase.cycles);
		add_graph_node(runtime, new_phase);

#ifdef NEW_DEBUG
		printf("\nadded new phase\n");
		print_profile(&new_phase->workload);
#endif
	}

#ifdef NEW_DEBUG
	printf("\ndumping all phases\n");
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		print_profile(&phases[i].workload);
	}
	printf("number of phases %d\n\n", runtime->classifier->numphases);
#endif

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
	if (profile->frq_target == runtime->sys->max_pstate)
	{
		// we can't duty cycle because we are at max freq
		set_perf(runtime, (uint16_t) profile->frq_target);
	}
	else
	{
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
}

struct phase_profile *update_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, struct phase_profile *prof, double phase_cycles)
{
	prof->phase_occurrences++;
	uint64_t prev_occurrences = prof->workload.occurrences;
	uint64_t new_occurrences = this_profile->occurrences;

	prof->workload.occurrences = prev_occurrences + new_occurrences;
	double prev_ratio = ((double) prev_occurrences / (double) prof->workload.occurrences);
	double new_ratio = ((double) new_occurrences / (double) prof->workload.occurrences);

	prof->workload.ipc = (prof->workload.ipc * prev_ratio) + (this_profile->ipc * new_ratio);
	prof->workload.mpc = (prof->workload.mpc * prev_ratio) + (this_profile->mpc * new_ratio);
	prof->workload.rpc = (prof->workload.rpc * prev_ratio) + (this_profile->rpc * new_ratio);
	prof->workload.epc = (prof->workload.epc * prev_ratio) + (this_profile->epc * new_ratio);
	prof->workload.bpc = (prof->workload.bpc * prev_ratio) + (this_profile->bpc * new_ratio);
	prof->workload.frq = (prof->workload.frq * prev_ratio) + (this_profile->frq * new_ratio);

	double prev_phase_ratio = ((double) prof->phase_occurrences - 1.0) /
			(double) prof->phase_occurrences;
	double new_phase_ratio = 1.0 / (double) prof->phase_occurrences;
	prof->cycles = prof->cycles * prev_phase_ratio + phase_cycles * new_phase_ratio;

	return prof;
}

struct phase_profile *add_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, double phase_cycles)
{
	struct phase_profile *newphase = 
		&runtime->classifier->phases[runtime->classifier->numphases];

	newphase->workload = *this_profile;
	newphase->cycles = phase_cycles;
	newphase->phase_occurrences = 1;

	runtime->classifier->numphases++;
	return newphase;
}


