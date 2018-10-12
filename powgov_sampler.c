#include <stdint.h>
#include <stdlib.h>
#include "powgov_sampler.h"
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"
#include "msr_counters.h"
#include "master.h"

int sample_data(struct powgov_runtime *runtime)
{
	// TODO: core stuff for multiple thread sampling
	int core = 0;
	if (core > runtime->cfg->threadcount || core < 0)
	{
		return -1;
	}
	
	uint64_t perf;
	uint64_t tsc;
	uint64_t energy;
	/* uint64_t rapl_throttled; */
	uint64_t therm;
	uint64_t perflimit;
	uint64_t instret;
	uint64_t aperf;
	uint64_t mperf;
	static char firstread = 1;
	struct pmc *pmcounters;
	pmc_storage(&pmcounters);
	read_msr_by_coord(0, core, 0, IA32_PERF_STATUS, &perf);
	read_msr_by_coord(0, core, 0, IA32_TIME_STAMP_COUNTER, &tsc);
	read_msr_by_coord(0, core, 0, MSR_IA32_APERF, &aperf);
	read_msr_by_coord(0, core, 0, MSR_IA32_MPERF, &mperf);
	read_msr_by_coord(0, core, 0, MSR_PKG_ENERGY_STATUS, &energy);
	//read_msr_by_coord(0, core, 0, MSR_PKG_PERF_STATUS, &rapl_throttled);
	read_msr_by_coord(0, core, 0, IA32_THERM_STATUS, &therm);
	read_msr_by_coord(0, core, 0, IA32_FIXED_CTR0, &instret);
	read_msr_by_coord(0, core, 0, MSR_CORE_PERF_LIMIT_REASONS, &perflimit);
	read_batch(COUNTERS_DATA);
	unsigned long idx = runtime->sampler->samplectrs[core];
	if (runtime->cfg->experimental)
	{
		runtime->sampler->l1->prev_sample = runtime->sampler->l1->new_sample;
		runtime->sampler->thread_samples[core][idx].frq_data = perf;
		runtime->sampler->thread_samples[core][idx].tsc_data = tsc;
		runtime->sampler->thread_samples[core][idx].aperf = aperf;
		runtime->sampler->thread_samples[core][idx].mperf = mperf;
		runtime->sampler->thread_samples[core][idx].energy_data = energy & 0xFFFFFFFF;
		/* runtime->sampler->thread_samples[core][idx].rapl_throttled = rapl_throttled & 0xFFFFFFFF; */
		runtime->sampler->thread_samples[core][idx].therm = therm;
		runtime->sampler->thread_samples[core][idx].perflimit = perflimit;
		runtime->sampler->thread_samples[core][idx].instret = instret;
		runtime->sampler->thread_samples[core][idx].llcmiss = *pmcounters->pmc0[0];
		runtime->sampler->thread_samples[core][idx].restalls = *pmcounters->pmc1[0];
		runtime->sampler->thread_samples[core][idx].exstalls = *pmcounters->pmc2[0];
		runtime->sampler->thread_samples[core][idx].branchret = *pmcounters->pmc3[0];
		runtime->sampler->l1->new_sample = runtime->sampler->thread_samples[core][idx];
		if (runtime->sampler->samplectrs[core] >= runtime->sampler->numsamples)
		{
			// every N samples do a buffered write to data file
			dump_data(runtime, runtime->files->sampler_dumpfiles);
			runtime->sampler->samplectrs[core] = 0;
		}
	}
	else
	{
		runtime->sampler->l1->prev_sample = runtime->sampler->l1->new_sample;
		runtime->sampler->l1->new_sample.frq_data = perf;
		runtime->sampler->l1->new_sample.tsc_data = tsc;
		runtime->sampler->l1->new_sample.energy_data = energy & 0xFFFFFFFF;
		/* runtime->sampler->l1->new_sample.rapl_throttled = rapl_throttled & 0xFFFFFFFF; */
		runtime->sampler->l1->new_sample.therm = therm;
		runtime->sampler->l1->new_sample.perflimit = perflimit;
		runtime->sampler->l1->new_sample.instret = instret;
		runtime->sampler->l1->new_sample.llcmiss = *pmcounters->pmc0[0];
		runtime->sampler->l1->new_sample.restalls = *pmcounters->pmc1[0];
		runtime->sampler->l1->new_sample.exstalls = *pmcounters->pmc2[0];
		runtime->sampler->l1->new_sample.branchret = *pmcounters->pmc3[0];
	}
	if (firstread)
	{
		runtime->sampler->first_sample = runtime->sampler->l1->new_sample;
		runtime->sampler->l2->new_sample = runtime->sampler->first_sample;
		runtime->sampler->l3->new_sample = runtime->sampler->first_sample;
		runtime->sampler->l1->prev_sample = runtime->sampler->l1->new_sample;
		/* runtime->sampler->l3->last_cyc = runtime->sampler->l1->new_sample.tsc_data; */
		firstread = 0;
	}
	if (runtime->sampler->total_samples % runtime->sampler->l2->interval == 0)
	{
		l2_analysis(runtime);
	}
	if (runtime->sampler->total_samples > 0 && 
			runtime->sampler->total_samples % runtime->sampler->l3->interval == 0)
	{
		l3_analysis(runtime);
	}

	runtime->sampler->samplectrs[core]++;
	runtime->sampler->total_samples++; // TODO: not thread safe
	return 0;
}

void init_sampling(struct powgov_runtime *runtime)
{
	runtime->sampler->thread_samples = (struct data_sample **) calloc(runtime->cfg->threadcount, sizeof(struct data_sample *));
	runtime->files->sampler_dumpfiles = (FILE **) calloc(runtime->cfg->threadcount, sizeof(FILE *));
	runtime->sampler->numsamples = runtime->sampler->sps;
	char fname[FNAMESIZE];
	int i;
	for (i = 0; i < runtime->cfg->threadcount; i++)
	{
		runtime->sampler->thread_samples[i] = (struct data_sample *) 
			calloc((runtime->sampler->numsamples) + 1, sizeof(struct data_sample));
		if (runtime->sampler->thread_samples[i] == NULL)
		{
			fprintf(stderr, "ERROR: out of memory\n");
			exit(-1);
		}
		snprintf((char *) fname, FNAMESIZE, "core%d.msrdat", i);
		runtime->files->sampler_dumpfiles[i] = fopen(fname, "w");
	}
}
