#include <math.h>
#include <string.h>
#include "powgov_profiles.h"
#include "powgov_util.h"
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"

double metric_distance(struct phase_profile *old, struct phase_profile *new, struct phase_profile *maximums,
	struct phase_profile *minimums)
{
	double ipcnorm = ((old->ipc - minimums->ipc) / (maximums->ipc - minimums->ipc)) -
		((new->ipc - minimums->ipc) / (maximums->ipc - minimums->ipc));

	double epcnorm = ((old->epc - minimums->epc) / (maximums->epc - minimums->epc)) -
		((new->epc - minimums->epc) / (maximums->epc - minimums->epc));
	//return sqrt(pow(ipcnorm, 2.0) + pow(mpcnorm, 2.0) + pow(rpcnorm, 2.0) + pow(epcnorm, 2.0) + pow(bpcnorm, 2.0));
	return sqrt(pow(ipcnorm, 2.0) + pow(epcnorm, 2.0));
}

// TODO: weighted averaging should scale first
void agglomerate_profiles(struct powgov_runtime *runtime)
{
	struct phase_profile old_profiles[MAX_PROFILES];
	int newidx = 0;
	struct phase_profile *profiles = runtime->classifier->profiles;
	memcpy(old_profiles, profiles, 
			runtime->classifier->numphases * sizeof(struct phase_profile));

	//double dist[MAX_PROFILES];
	//memset(dist, 0, sizeof(double) * MAX_PROFILES);

	char valid[MAX_PROFILES];
	memset(valid, 1, MAX_PROFILES);

	int numcombines = 0;

	//printf("(glom) num phases %d\n", runtime->classifier->numphases);

	int i;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		if (valid[i] == 0)
		{
			continue;
		}

		char matches[MAX_PROFILES];
		memset(matches, 0, MAX_PROFILES);
		uint64_t occurrence_sum = old_profiles[i].occurrences;
		double cycles_sum = old_profiles[i].avg_cycle;

		int j;
		for (j = 0; j < runtime->classifier->numphases; j++)
		{
			if (i != j && valid[j])
			{
				double dist = metric_distance(&old_profiles[i], &old_profiles[j], 
						&runtime->classifier->prof_maximums, &runtime->classifier->prof_minimums);
				if (dist < runtime->classifier->dist_thresh)
				{
					matches[j] = 1;
					matches[i] = 1;
					valid[j] = 0;
					valid[i] = 0;
					occurrence_sum += old_profiles[j].occurrences;
					cycles_sum += old_profiles[j].occurrences;
					numcombines++;
				}
			}
		}
		if (matches[i])
		{
			memcpy(&runtime->classifier->profiles[newidx], &old_profiles[i], sizeof(struct phase_profile));
			profiles[newidx].ipc *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].mpc *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].rpc *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].epc *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].bpc *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].num_throttles *=
					(profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].avg_cycle = cycles_sum / numcombines;
			//profiles[newidx].avg_cycle *= (profiles[newidx].occurrences / (double) occurrence_sum);
			profiles[newidx].avg_frq *= (profiles[newidx].occurrences / (double) occurrence_sum);
			int numglom = 1;
			for (j = 0; j < runtime->classifier->numphases; j++)
			{
				if (i != j && matches[j])
				{
					struct phase_profile scaled_profile;
					frequency_scale_phase(&old_profiles[j], old_profiles[j].avg_frq,
							old_profiles[i].avg_frq, &scaled_profile);
					//printf("(glom) combining %d and %d\n", i, j);
					if (j == runtime->classifier->recentphase)
					{
						runtime->classifier->recentphase = newidx;
					}
					profiles[newidx].ipc += scaled_profile.ipc *
							(scaled_profile.occurrences / (double) occurrence_sum);
					profiles[newidx].mpc += scaled_profile.mpc *
							(scaled_profile.occurrences / (double) occurrence_sum);
					profiles[newidx].rpc += scaled_profile.rpc *
							(scaled_profile.occurrences / (double) occurrence_sum);
					profiles[newidx].epc += scaled_profile.epc *
							(scaled_profile.occurrences / (double) occurrence_sum);
					profiles[newidx].bpc += scaled_profile.bpc *
							(scaled_profile.occurrences / (double) occurrence_sum);
					profiles[newidx].num_throttles += scaled_profile.num_throttles *
							(scaled_profile.occurrences / (double) occurrence_sum);
					/* profiles[newidx].avg_cycle += scaled_profile.avg_cycle * */
					/* 		(scaled_profile.occurrences / (double) occurrence_sum); */

					if (scaled_profile.frq_high > profiles[newidx].frq_high)
					{
						profiles[newidx].frq_high = scaled_profile.frq_high;
					}
					if (scaled_profile.frq_low < profiles[newidx].frq_low)
					{
						profiles[newidx].frq_low = scaled_profile.frq_low;
					}
					profiles[newidx].occurrences += scaled_profile.occurrences;
					profiles[newidx].avg_frq += scaled_profile.avg_frq * (scaled_profile.occurrences / (double) occurrence_sum);
					numglom++;
				}
			}
			// TODO: instead of clearing, update prev_phases
			profiles[newidx].lastprev = 0;
			/* memset(profiles[newidx].prev_phases, -1, MAX_HISTORY); */
			newidx++;
		}
	}
	if (numcombines > 0)
	{
		// slide everything over that wasn't combined
		for (i = 0; i < runtime->classifier->numphases; i++)
		{
			if (valid[i] && i == newidx)
			{
				//printf("(glom) skipping id %d\n", i);
				newidx++;
			}
			else if (valid[i])
			{
				//printf("(glom) sliding %d to %d\n", i, newidx);
				if (i == runtime->classifier->recentphase)
				{
					runtime->classifier->recentphase = newidx;
				}
				memcpy(&profiles[newidx], &old_profiles[i], sizeof(struct phase_profile));
				runtime->classifier->profiles[newidx].lastprev = 0;
				/* memset(profiles[newidx].prev_phases, -1, MAX_HISTORY); */

				newidx++;
			}
		}
		runtime->classifier->numphases = newidx;
	}
	//printf("(glom) runtime->classifier->recentphase is now %d\n", runtime->classifier->recentphase);
}

void remove_unused(struct powgov_runtime *runtime)
{
	if (runtime->classifier->numphases <= 0)
	{
		return;
	}
	//printf("(remove) numphases %d\n", numphases);

	int i;
	//uint64_t occurrence_sum = 0;
	char valid[MAX_PROFILES];
	memset(valid, 1, MAX_PROFILES);
	int numinvalid = 0;
	int firstinvalid = -1;
	struct phase_profile *profiles = runtime->classifier->profiles;

	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		//occurrence_sum += old_profiles[i].occurrences;
		if (profiles[i].occurrences <= 1)
		{
			//printf("(remove) removing id %d with %lu occurrences\n", i, profiles[i].occurrences);
			valid[i] = 0;
			numinvalid++;
			if (firstinvalid < 0)
			{
				firstinvalid = i;
			}
		}
	}
	/*
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		if (old_profiles[i].occurrences / (double) occurrence_sum < PRUNE_THRESH)
		{
			valid[i] = 0;
			numinvalid++;
		}
	}
	*/
	if (numinvalid == 0)
	{
		return;
	}

	struct phase_profile old_profiles[MAX_PROFILES];
	int newidx = 0;
	memcpy(old_profiles, profiles, runtime->classifier->numphases * sizeof(struct phase_profile));

	newidx = firstinvalid;
	for (i = firstinvalid; i < runtime->classifier->numphases; i++)
	{
		if (i == newidx && valid[i])
		{
			//printf("(remove) skipping id %d\n", i);
			newidx++;
		}
		else if (valid[i])
		{
			//printf("(remove) sliding id %d to %d\n", i, newidx);
			if (i == runtime->classifier->recentphase)
			{
				runtime->classifier->recentphase = newidx;
			}
			memcpy(&profiles[newidx], &old_profiles[i], sizeof(struct phase_profile));
			newidx++;
		}
	}
 runtime->classifier->numphases = newidx;
	//printf("(remove) runtime->classifier->numphases is now %d\n", runtime->classifier->numphases);
	//printf("(remove) runtime->classifier->recentphase is now %d\n", runtime->classifier->recentphase);
}

void update_minmax(struct phase_profile *this_profile, struct phase_profile *maximums, 
		struct phase_profile *minimums)
{
	if (this_profile->ipc > maximums->ipc)
	{
		maximums->ipc = this_profile->ipc;
	}
	if (this_profile->mpc > maximums->mpc)
	{
		maximums->mpc = this_profile->mpc;
	}
	if (this_profile->rpc > maximums->rpc)
	{
		maximums->rpc = this_profile->rpc;
	}
	if (this_profile->epc > maximums->epc)
	{
		maximums->epc = this_profile->epc;
	}
	if (this_profile->bpc > maximums->bpc)
	{
		maximums->bpc = this_profile->bpc;
	}

	if (this_profile->ipc < minimums->ipc)
	{
		minimums->ipc = this_profile->ipc;
	}
	if (this_profile->mpc < minimums->mpc)
	{
		minimums->mpc = this_profile->mpc;
	}
	if (this_profile->rpc < minimums->rpc)
	{
		minimums->rpc = this_profile->rpc;
	}
	if (this_profile->epc < minimums->epc)
	{
		minimums->epc = this_profile->epc;
	}
	if (this_profile->bpc < minimums->bpc)
	{
		minimums->bpc = this_profile->bpc;
	}
}

void print_profile(struct phase_profile *prof)
{
	printf("ipc: %lf\nmpc %lf\nrpc %lf\nepc %lf\nbpc %lf\n", prof->ipc, prof->mpc,
			prof->rpc, prof->epc, prof->bpc);
}

void update_profile(struct powgov_runtime *runtime, struct phase_profile *this_profile, int profidx, uint64_t perf, unsigned this_throttle, double avgfrq, int lastphase)
{
	if (profidx > runtime->classifier->numphases)
	{
		printf("ERROR: profile does not exist\n");
		return;
	}
	struct phase_profile *profiles = runtime->classifier->profiles;

	profiles[profidx].ipc = (profiles[profidx].ipc *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(this_profile->ipc * (1.0 / (profiles[profidx].occurrences + 1.0)));

	profiles[profidx].mpc = (profiles[profidx].mpc *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(this_profile->mpc * (1.0 / (profiles[profidx].occurrences + 1.0)));

	profiles[profidx].rpc = (profiles[profidx].rpc *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(this_profile->rpc * (1.0 / (profiles[profidx].occurrences + 1.0)));

	profiles[profidx].epc = (profiles[profidx].epc *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(this_profile->epc * (1.0 / (profiles[profidx].occurrences + 1.0)));

	profiles[profidx].bpc = (profiles[profidx].bpc *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(this_profile->bpc * (1.0 / (profiles[profidx].occurrences + 1.0)));

//	profiles[profidx].num_throttles = (profiles[profidx].num_throttles *
//			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
//			(this_profile->num_throttles * (1.0 / (profiles[profidx].occurrences + 1.0)));
	profiles[profidx].num_throttles += this_throttle;

	if (perf > profiles[profidx].frq_high)
	{
		profiles[profidx].frq_high = perf;
	}
	if (perf < profiles[profidx].frq_low)
	{
		profiles[profidx].frq_low = perf;
	}

	profiles[profidx].avg_frq = (profiles[profidx].avg_frq *
			(profiles[profidx].occurrences / (profiles[profidx].occurrences + 1.0))) +
			(avgfrq * (1.0 / (profiles[profidx].occurrences + 1.0)));
	profiles[profidx].occurrences++;

	/* char last = profiles[profidx].lastprev; */
	/* int i; */
	/* for (i = 0; i <= last && i < MAX_HISTORY; i++) */
	/* { */
	/* 	if (profiles[profidx].prev_phases[i] == lastphase) */
	/* 	{ */
	/* 		// it's already there */
	/* 		return; */
	/* 	} */
	/* } */
	/* if (last < MAX_HISTORY) */
	/* { */
	/* 	last++; */
	/* 	profiles[profidx].prev_phases[last] = lastphase; */
	/* 	profiles[profidx].lastprev = last; */
	/* } */
}

void add_profile(struct powgov_runtime *runtime, struct phase_profile *this_profile, uint64_t perf, unsigned this_throttle, double avgfrq, int lastphase)
{
	int numphases = runtime->classifier->numphases;
	struct phase_profile *profiles = runtime->classifier->profiles;
	profiles[numphases].ipc = this_profile->ipc;
	profiles[numphases].mpc = this_profile->mpc;
	profiles[numphases].rpc = this_profile->rpc;
	profiles[numphases].epc = this_profile->epc;
	profiles[numphases].bpc = this_profile->bpc;

	profiles[numphases].avg_frq = avgfrq;
	profiles[numphases].frq_high = runtime->sys->min_pstate;
	profiles[numphases].frq_low = runtime->sys->max_pstate;
	profiles[numphases].frq_target = (double) runtime->sys->max_pstate;
	profiles[numphases].avg_cycle = 0;
	profiles[numphases].num_throttles = this_throttle;
	profiles[numphases].occurrences = 0;
	/* memset(profiles[numphases].prev_phases, -1, MAX_HISTORY); */
	/* profiles[numphases].prev_phases[0] = lastphase; */
	profiles[numphases].lastprev = 0;
	profiles[numphases].class = 4;
	profiles[numphases].unthrot_count = 0;
	profiles[numphases].reclass_count = RECLASSIFY_INTERVAL;
	profiles[numphases].frq_duty_count = 0;
	/* profiles[numphases].mem_fn1_ctr = 0; */
	/* profiles[numphases].mem_fn2_ctr = 0; */
	/* profiles[numphases].mem_fn3_ctr = 0; */

	 runtime->classifier->numphases++;
#ifdef DEBUG
	printf("Added new phase profile %d\n", runtime->classifier->numphases);
#endif
}

int classify_phase(struct powgov_runtime *runtime, struct phase_profile *phase, uint64_t perf)
{
	int i = -1;
	int minidx = -1;
	double mindist = DBL_MAX;
	double freq = ((double) perf) / 10.0;
	struct phase_profile *prof_class = runtime->classifier->prof_class;
	prof_class[0].ipc = freq * CLASS_CPU_SLOPE_IPC + CLASS_CPU_INTERCEPT_IPC;
	prof_class[0].epc = freq * CLASS_CPU_SLOPE_EPC + CLASS_CPU_INTERCEPT_EPC;
	prof_class[1].ipc = freq * CLASS_MEM_SLOPE_IPC + CLASS_MEM_INTERCEPT_IPC;
	prof_class[1].epc = freq * CLASS_MEM_SLOPE_EPC + CLASS_MEM_INTERCEPT_EPC;
	// skip IO class, it doesn't change much
	prof_class[3].ipc = freq * CLASS_MIX_SLOPE_IPC + CLASS_MIX_INTERCEPT_IPC;
	prof_class[3].epc = freq * CLASS_MIX_SLOPE_EPC + CLASS_MIX_INTERCEPT_EPC;
	for (i = 0; i < NUM_CLASSES; i++)
	{
		double dist = metric_distance(phase, &prof_class[i], &runtime->classifier->prof_maximums, &runtime->classifier->prof_minimums);
		if (dist < mindist)
		{
			mindist = dist;
			minidx = i;
		}
	}
	phase->class = minidx;
	phase->reclass_count = 0;
	if (runtime->sampler->l3->seq_end < 0)
	{
		runtime->sampler->l3->seq_end = 0;
		runtime->sampler->l3->sequence[runtime->sampler->l3->seq_end] = phase->class;
		runtime->sampler->l3->seq_cycles[runtime->sampler->l3->seq_end] =
				runtime->sampler->l1->new_sample.tsc_data - 
				runtime->sampler->l3->last_cyc;
		runtime->sampler->l3->last_cyc = runtime->sampler->l1->new_sample.tsc_data;
	}
	else if (phase->class != runtime->sampler->l3->sequence[runtime->sampler->l3->seq_end])
	{
		if (runtime->sampler->l3->seq_end < MAX_L3_SEQ)
		{
			runtime->sampler->l3->seq_end++;
			runtime->sampler->l3->sequence[runtime->sampler->l3->seq_end] = phase->class;
			runtime->sampler->l3->seq_cycles[runtime->sampler->l3->seq_end] =
					runtime->sampler->l1->new_sample.tsc_data - 
					runtime->sampler->l3->last_cyc;
			runtime->sampler->l3->last_cyc = runtime->sampler->l1->new_sample.tsc_data;
		}
	}
#ifdef DEBUG
	if (phase->epc > 2.5)
	{
		printf("CLASS %d freq %lf\n", phase->class, freq);
		printf("phase ipc %lf epc %lf\n", phase->ipc, phase->epc);
		printf("cpu scale ipc %lf epc %lf\n", prof_class[0].ipc, prof_class[0].epc);
		printf("mem scale ipc %lf epc %lf\n", prof_class[1].ipc, prof_class[1].epc);
		printf("mix scale ipc %lf epc %lf\n", prof_class[2].ipc, prof_class[2].epc);
	}
#endif
	return minidx;
}

void classify_and_react(struct powgov_runtime *runtime, int phase, char wasthrottled, uint64_t perf)
{
	int class;
	struct phase_profile *profiles = runtime->classifier->profiles;
	// avoid reclassifying every timestep
	if (profiles[phase].reclass_count >= RECLASSIFY_INTERVAL)
	{
		class = classify_phase(runtime, &profiles[phase], perf);
	}
	else
	{
		class = profiles[phase].class;
		profiles[phase].reclass_count++;
	}
	if (runtime->cfg->man_cpu_ctrl)
	{
		if (class == CLASS_CPU)
		{
			if (runtime->cfg->cpu_frq_override != 0.0)
			{
				profiles[phase].frq_target = runtime->cfg->cpu_frq_override * 10.0;
				set_perf(runtime, (uint16_t) profiles[phase].frq_target);
				return;
			}
			if (runtime->power->excursion)
			{
				profiles[phase].frq_target -= runtime->cfg->frq_change_step * 2;
				runtime->power->excursion = 0;
			}
		}
		if (wasthrottled && runtime->cfg->throttle_avoid)
		{
			set_perf(runtime, (uint16_t) (profiles[phase].frq_target - runtime->cfg->frq_change_step * 
						runtime->cfg->throttle_step_mult));
			if (perf < profiles[phase].frq_target)
			{
				profiles[phase].frq_target -= runtime->cfg->frq_change_step * runtime->cfg->throttle_step_mult;
			}
			profiles[phase].unthrot_count = 0;
		}
		else if (class == CLASS_MEM)
		{
			if (runtime->cfg->mem_pow_shift)
			{
				profiles[phase].frq_target = runtime->classifier->mem_freq_throttle;
				set_perf(runtime, (uint16_t) profiles[phase].frq_target);
				// power shifting heuristic
				/* if (mem_fn1_ctr < MEM_FN1_LIM) */
				/* { */
				/* 	profiles[phase].frq_target = (double) mem_throttle_low; */
				/* 	mem_fn1_ctr++; */
				/* 	set_perf(runtime, (uint16_t) profiles[phase].frq_target); */
				/* } */
				/* else if (mem_fn2_ctr < MEM_FN2_LIM) */
				/* { */
				/* 	profiles[phase].frq_target = (double)  */
				/* 		(mem_throttle_low + MEM_STEP1 <= runtime->sys->max_pstate ? */
				/* 		 mem_throttle_low + MEM_STEP1 : runtime->sys->max_pstate); */
				/* 	mem_fn2_ctr++; */
				/* 	set_perf(runtime, (uint16_t) profiles[phase].frq_target); */
				/* } */
				/* else if (mem_fn3_ctr < MEM_FN3_LIM) */
				/* { */
				/* 	profiles[phase].frq_target = (double) */
				/* 		(mem_throttle_low + MEM_STEP2 <= runtime->sys->max_pstate ? */
				/* 		 mem_throttle_low + MEM_STEP2 : runtime->sys->max_pstate); */
				/* 	mem_fn3_ctr++; */
				/* 	set_perf(runtime, (uint16_t) profiles[phase].frq_target); */
				/* } */
				/* else */
				/* { */
				/* 	profiles[phase].frq_target = (double) runtime->sys->max_pstate; */
				/* 	set_perf(runtime, (uint16_t) profiles[phase].frq_target); */
				/* } */
			}
			else if (runtime->cfg->mem_frq_override != 0.0)
			{
				// power shifting with user override
				profiles[phase].frq_target = runtime->cfg->mem_frq_override * 10.0;
				set_perf(runtime, (uint16_t) profiles[phase].frq_target);
			}
			else
			{
				// default behavior
				profiles[phase].frq_target = (double) runtime->sys->max_pstate;
				set_perf(runtime, (uint16_t) profiles[phase].frq_target);
			}
		}
		else
		{
			// TODO fix this control flow nightmare
			if (perf > profiles[phase].frq_target && profiles[phase].frq_target < runtime->sys->max_pstate)
			{
				profiles[phase].frq_target += runtime->cfg->frq_change_step;
			}
			if (profiles[phase].frq_target < runtime->sys->max_pstate &&
				profiles[phase].unthrot_count >= runtime->cfg->fup_timeout)
			{
				profiles[phase].frq_target += runtime->cfg->frq_change_step;
				// TODO may want to have separate counter for this, timeout before each freq increase
				// although not wrong because self throttling...
				profiles[phase].unthrot_count = 0;
			}
			else
			{
				profiles[phase].unthrot_count++;
			}
			double ratio = profiles[phase].frq_target - (double)((uint16_t) profiles[phase].frq_target);
			profiles[phase].frq_duty_count++;
			if (profiles[phase].frq_duty_count < ratio * runtime->cfg->frq_duty_length)
			{
				// set perf to ciel
				set_perf(runtime, (uint16_t) profiles[phase].frq_target + 1);
			}
			else
			{
				// set perf to floor
				set_perf(runtime, (uint16_t) profiles[phase].frq_target);
			}
			if (profiles[phase].frq_duty_count > runtime->cfg->frq_duty_length)
			{
				profiles[phase].frq_duty_count = 0;
			}
			/* set_perf(runtime, profiles[phase].frq_target); */
		}
	}
}

// TODO: unused func
double ipc_scale(double ipc_unscaled, double frq_source, double frq_target)
{
	if (frq_target + SCALE_THRESH > frq_source && frq_target - SCALE_THRESH < frq_source)
	{
		return ipc_unscaled;
	}
	double cpuipc = frq_source * CLASS_CPU_SLOPE_IPC + CLASS_CPU_INTERCEPT_IPC;
	double memipc = frq_source * CLASS_MEM_SLOPE_IPC + CLASS_MEM_INTERCEPT_IPC;
	double ipc_percent = (ipc_unscaled - memipc) / (cpuipc - memipc);
	double result = ((ipc_percent * CLASS_CPU_SLOPE_IPC) + ((1.0 - ipc_percent) * CLASS_MEM_SLOPE_IPC)) *
		frq_target + (ipc_percent * CLASS_CPU_INTERCEPT_IPC) + ((1.0 - ipc_percent) * CLASS_MEM_SLOPE_IPC);

	if (result < 0.0)
	{
		result = 0.0;
	}

	return result;
}

// TODO: keep working on scaling accuracy
void frequency_scale_phase(struct phase_profile *unscaled_profile, double frq_source, double frq_target, struct phase_profile *scaled_profile)
{
	*scaled_profile = *unscaled_profile;
	// if the frequencies are already close then just return the copy, don't scale
	if (frq_target + SCALE_THRESH > frq_source && frq_target - SCALE_THRESH < frq_source)
	{
		return;
	}
	double cpuipc = frq_source * CLASS_CPU_SLOPE_IPC + CLASS_CPU_INTERCEPT_IPC;
	double memipc = frq_source * CLASS_MEM_SLOPE_IPC + CLASS_MEM_INTERCEPT_IPC;
	double ipc_percent = (unscaled_profile->ipc - memipc) / (cpuipc - memipc);

	double cpuepc = frq_source * CLASS_CPU_SLOPE_EPC + CLASS_CPU_INTERCEPT_EPC;
	double memepc = frq_source * CLASS_MEM_SLOPE_EPC + CLASS_MEM_INTERCEPT_EPC;
	// cpu and mem switched here because CPU has minimum in this case
	double epc_percent = (unscaled_profile->epc - cpuepc) / (memepc - cpuepc);

	// if the percentages do not add close to 1, then this data point is atypical so don't scale it
	if (ipc_percent + epc_percent < SCALE_OUTLIER_THRESH_LOW || ipc_percent + epc_percent > SCALE_OUTLIER_THRESH_HIGH)
	{
		return;
	}
	scaled_profile->ipc = ((ipc_percent * CLASS_CPU_SLOPE_IPC) + ((1.0 - ipc_percent) * CLASS_MEM_SLOPE_IPC)) *
		frq_target + (ipc_percent * CLASS_CPU_INTERCEPT_IPC) + ((1.0 - ipc_percent) * CLASS_MEM_SLOPE_IPC);
	scaled_profile->epc = ((epc_percent * CLASS_MEM_SLOPE_EPC) + ((1.0 - epc_percent) * CLASS_CPU_SLOPE_EPC)) *
		frq_target + (epc_percent * CLASS_MEM_INTERCEPT_EPC) + ((1.0 - epc_percent) * CLASS_CPU_SLOPE_EPC);
	if (scaled_profile->ipc < 0.0)
	{
		scaled_profile->ipc = 0.0;
	}
	if (scaled_profile->epc < 0.0)
	{
		scaled_profile->epc = 0.0;
	}
#ifdef DEBUG
	printf("\nfrq source %lf, frq target %lf\n", frq_source, frq_target);
	printf("ipc scale from %lf to %lf\nepc scaled from %lf to %lf\n", unscaled_profile->ipc, scaled_profile->ipc,
			unscaled_profile->epc, scaled_profile->epc);
	printf("ipc percent %lf, epc percent %lf\n", ipc_percent, epc_percent);
#endif
}
