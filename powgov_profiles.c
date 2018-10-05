#include <math.h>
#include <float.h>
#include <string.h>
#include "powgov_profiles.h"
#include "powgov_util.h"
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"
// TODO: occurrences does not mean the same thing anymore, re-work weighted averages

double workload_metric_distance(struct workload_profile *old, struct workload_profile *new, struct workload_profile *maximums, struct workload_profile *minimums)
{
	double ipcnorm = ((old->ipc - minimums->ipc) / (maximums->ipc - minimums->ipc)) -
		((new->ipc - minimums->ipc) / (maximums->ipc - minimums->ipc));

	double epcnorm = ((old->epc - minimums->epc) / (maximums->epc - minimums->epc)) -
		((new->epc - minimums->epc) / (maximums->epc - minimums->epc));
	//return sqrt(pow(ipcnorm, 2.0) + pow(mpcnorm, 2.0) + pow(rpcnorm, 2.0) + pow(epcnorm, 2.0) + pow(bpcnorm, 2.0));
	return sqrt(pow(ipcnorm, 2.0) + pow(epcnorm, 2.0));
}

double phase_metric_distance(struct phase_profile *old, struct phase_profile *new, struct workload_profile *maximums, struct workload_profile *minimums)
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
	struct workload_profile old_profiles[MAX_PROFILES];
	int newidx = 0;
	struct workload_profile *profiles = runtime->classifier->profiles;
	memcpy(old_profiles, profiles, 
			runtime->classifier->numphases * sizeof(struct workload_profile));

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
			memcpy(&runtime->classifier->profiles[newidx], &old_profiles[i], sizeof(struct workload_profile));
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
					struct workload_profile scaled_profile;
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
				memcpy(&profiles[newidx], &old_profiles[i], sizeof(struct workload_profile));
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
	struct workload_profile *profiles = runtime->classifier->profiles;

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

	struct workload_profile old_profiles[MAX_PROFILES];
	int newidx = 0;
	memcpy(old_profiles, profiles, runtime->classifier->numphases * sizeof(struct workload_profile));

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
			memcpy(&profiles[newidx], &old_profiles[i], sizeof(struct workload_profile));
			newidx++;
		}
	}
 runtime->classifier->numphases = newidx;
	//printf("(remove) runtime->classifier->numphases is now %d\n", runtime->classifier->numphases);
	//printf("(remove) runtime->classifier->recentphase is now %d\n", runtime->classifier->recentphase);
}

void update_minmax(struct powgov_runtime *runtime, struct workload_profile *this_profile)
{
	if (this_profile->ipc > runtime->classifier->prof_maximums.ipc)
	{
		runtime->classifier->prof_maximums.ipc = this_profile->ipc;
	}
	if (this_profile->mpc > runtime->classifier->prof_maximums.mpc)
	{
		runtime->classifier->prof_maximums.mpc = this_profile->mpc;
	}
	if (this_profile->rpc > runtime->classifier->prof_maximums.rpc)
	{
		runtime->classifier->prof_maximums.rpc = this_profile->rpc;
	}
	if (this_profile->epc > runtime->classifier->prof_maximums.epc)
	{
		runtime->classifier->prof_maximums.epc = this_profile->epc;
	}
	if (this_profile->bpc > runtime->classifier->prof_maximums.bpc)
	{
		runtime->classifier->prof_maximums.bpc = this_profile->bpc;
	}

	if (this_profile->ipc < runtime->classifier->prof_minimums.ipc)
	{
		runtime->classifier->prof_minimums.ipc = this_profile->ipc;
	}
	if (this_profile->mpc < runtime->classifier->prof_minimums.mpc)
	{
		runtime->classifier->prof_minimums.mpc = this_profile->mpc;
	}
	if (this_profile->rpc < runtime->classifier->prof_minimums.rpc)
	{
		runtime->classifier->prof_minimums.rpc = this_profile->rpc;
	}
	if (this_profile->epc < runtime->classifier->prof_minimums.epc)
	{
		runtime->classifier->prof_minimums.epc = this_profile->epc;
	}
	if (this_profile->bpc < runtime->classifier->prof_minimums.bpc)
	{
		runtime->classifier->prof_minimums.bpc = this_profile->bpc;
	}
}

void print_profile(struct workload_profile *prof)
{
	printf("ipc: %lf\nmpc %lf\nrpc %lf\nepc %lf\nbpc %lf\n", prof->ipc, prof->mpc,
			prof->rpc, prof->epc, prof->bpc);
}

void update_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, struct phase_profile *prof, double phase_cycles)
{
	prof->phase_occurrences++;
	uint64_t prev_occurrences = prof->workload.occurrences;
	uint64_t new_occurrences = this_profile->occurrences;

	prof->workload.occurrences = prev_occurrences + new_occurrences;
	prof->workload.num_throttles += this_profile->num_throttles;

	prof->workload.ipc = (prof->workload.ipc *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->ipc * (new_occurrences / prof->workload.occurrences));
	prof->workload.mpc = (prof->workload.mpc *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->mpc * (new_occurrences / prof->workload.occurrences));
	prof->workload.rpc = (prof->workload.rpc *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->rpc * (new_occurrences / prof->workload.occurrences));
	prof->workload.epc = (prof->workload.epc *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->epc * (new_occurrences / prof->workload.occurrences));
	prof->workload.bpc = (prof->workload.bpc *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->bpc * (new_occurrences / prof->workload.occurrences));
	prof->workload.frq = (prof->workload.frq *
			(prev_occurrences / prof->workload.occurrences)) +
			(this_profile->frq * (new_occurrences / prof->workload.occurrences));

	prof->cycles = (prof->cycles * (prof->phase_occurrences - 1.0) +
			phase_cycles * (1.0 / prof->phase_occurrences));
}

void add_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, double phase_cycles)
{
	struct phase_profile *newphase = 
		runtime->classifier->phases[runtime->classifier->numphases];

	newphase->workload = this_profile;
	newphase->cycles = phase_cycles;
	newphase->phase_occurrences = 1;

	runtime->classifier->numphases++;
}

int classify_workload(struct powgov_runtime *runtime, struct workload_profile *phase, uint64_t perf)
{
	int i = -1;
	int minidx = -1;
	double mindist = DBL_MAX;
	double freq = ((double) perf) / 10.0;
	struct workload_profile *prof_class = runtime->classifier->prof_class;
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

// TODO: keep working on scaling accuracy
void frequency_scale_phase(struct workload_profile *unscaled_profile, double frq_source, double frq_target, struct workload_profile *scaled_profile)
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
