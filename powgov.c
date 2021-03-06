#define _GNU_SOURCE // required for sched
#define _BSD_SOURCE // required for usleep
#include <unistd.h>
#include <stdio.h>
#include <sched.h>
#include <signal.h>
#include <float.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "powgov.h"
#include "powgov_util.h"
#include "powgov_sampler.h"
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"
#include "cpuid.h"

#define HEADERSTRING "########################################"
#define FNAMESIZE 32
#define VERSIONSTRING "0.3"
#define ARG_ERROR {\
	if (j + 1 > argc - 1)\
	{\
		printf("Error: value required for argument %s\n", argv[j]);\
		exit(-1);\
	}\
}

// Globals
int LOOP_CTRL = 1;
double VERBOSE;
double REPORT;
char CLASS_NAMES[NUM_CLASSES + 1][8] = {"CPU\0", "MEM\0", "IO\0", "MIX\0", "UNK\0"};

// this can be used to dump a report to a file if the program has an error
void dump_error_report(struct powgov_runtime *runtime)
{
	char fname[32];
	snprintf(fname, 32, "powgov_error");
	FILE *error = fopen(fname, "w");
	double avgrate = 1.0; // does not mean anything for error reporting
	dump_phaseinfo(runtime, error, &avgrate);
	fclose(error);
}

// TODO: why is avgrate a pointer?
// this dumps the phase descriptors, and graph information
void dump_phaseinfo(struct powgov_runtime *runtime, FILE *outfile, double *avgrate)
{
	int i;
	uint64_t recorded_steps = 0;
	struct phase_profile *phases = runtime->classifier->phases;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		recorded_steps += phases[i].workload.occurrences;
	}
	double totaltime = 0.0;
	double totalpct = 0.0;
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		double pct = (double) phases[i].workload.occurrences / (double) recorded_steps;
		if (avgrate != NULL)
		{
			fprintf(outfile, "PHASE ID %4d\t %12.3lf seconds\t(%6.2lf%%)\t%16.0lf cycles\n",
					i, *avgrate * phases[i].workload.occurrences, pct * 100.0, phases[i].cycles);
			totaltime += *avgrate * phases[i].workload.occurrences;
			//totalpct += pct * 100.0;
		}
	}
	if (avgrate != NULL)
	{
		totalpct = (double) recorded_steps / (double) runtime->sampler->samplectrs[0];
		fprintf(outfile, "TOTAL\t\t%12.2lf\t(%5.2lf%% accounted)\n", totaltime, totalpct * 100.0);
	}
	
	for (i = 0; i < runtime->classifier->numphases; i++)
	{
		// ignore phases that are less than x% of program
		double pct = (double) phases[i].workload.occurrences / (double) recorded_steps;
		if (pct < runtime->classifier->pct_thresh)
		{
			continue;
		}
		if (avgrate != NULL)
		{
			fprintf(outfile, "\nPHASE ID %4d\t %12.3lf seconds\t(%6.2lf%%)\t%16.0lf cycles\n",
					i, *avgrate * phases[i].workload.occurrences, pct * 100.0, phases[i].cycles);
		}
		else
		{
			fprintf(outfile, "\nPHASE ID %d\t(%3.0lf%%)\n", i, pct * 100);
		}
		fprintf(outfile, "\tinstructions per cycle        %6.4lf\n", phases[i].workload.ipc);
		fprintf(outfile, "\tLLC misses per cycle          %6.4lf\n", phases[i].workload.mpc);
		fprintf(outfile, "\tresource stalls per cycle     %6.4lf\n", phases[i].workload.rpc);
		fprintf(outfile, "\texecution stalls per cycle    %6.4lf\n", phases[i].workload.epc);
		fprintf(outfile, "\tbranch instructions per cycle %6.4lf\n", phases[i].workload.bpc);
		fprintf(outfile, "\tphase occurrences             %lu\n", 
				phases[i].workload.occurrences);
		fprintf(outfile, "\tfrequency                     %1.2lf\n",
				FRQ_AS_GHZ(phases[i].workload.frq));
		fprintf(outfile, "\tfrq target                    %1.2lf\n",
				FRQ_AS_GHZ(phases[i].workload.frq_target));
		fprintf(outfile, "\tclass                         %s (#%d)\n",
				CLASS_NAMES[(int)phases[i].workload.class], phases[i].workload.class);

		struct l3_graph_node *itr;
		struct l3_graph_node *begin = runtime->sampler->l3->graph;
		for (itr = begin; itr < begin + MAX_L3_GRAPH && itr->phase != &phases[i]; itr++);
		fprintf(outfile, "\tpredecessors                  %u\n", itr->reference_counter);
		fprintf(outfile, "\tSuccessors:\n");
		int j;
		for (j = 0; j < 3; j++)
		{
			if (itr->next[j] != NULL)
			{
				fprintf(outfile, "\t\t%4d (occurs %16lu)\n",
						(int) (itr->next[j]->phase - runtime->classifier->phases),
						itr->next_occurrences[j]);
			}
		}
		fprintf(outfile, "\tFrequency Throttles:\n");
		for (j = 0; j < MAX_PSTATE_NUM; j++)
		{
			if (phases[i].workload.throttlehist[j] != 0)
			{
				fprintf(outfile, "\t\t%2.2f: %u\n", (j + 8) / 10.0,
						phases[i].workload.throttlehist[j]);
			}
		}

		fprintf(outfile, "\tSimilarity:\n");
		for (j = 0; j < runtime->classifier->numphases; j++)
		{
			double lpct = (double) phases[j].workload.occurrences / (double) recorded_steps;
			if (lpct > runtime->classifier->pct_thresh)
			{
				struct workload_profile scaled_profile;
				frequency_scale_phase(&phases[j].workload, phases[j].workload.frq, phases[i].workload.frq, &scaled_profile);
				double dist = workload_metric_distance(&scaled_profile, &phases[i].workload, &runtime->classifier->prof_maximums);
				fprintf(outfile, "\t\tdistance from %4d: %4.4lf (%16.0lf cycles)\n", j, dist,
						fabs(phases[i].cycles - phases[j].cycles));
			}
		}
	}
}

// this dumps the non-phase based data collected during the run
void dump_data(struct powgov_runtime *runtime, FILE **outfile)
{
	struct data_sample **thread_samples = runtime->sampler->thread_samples;
	static long last_dump_loc = 0;
	int j;
	for (j = 0; j < runtime->cfg->threadcount; j++)
	{
		if (last_dump_loc == 0)
		{
			fprintf(outfile[j], "freq\tp-state\ttsc\tpower\trapl-throttle-cycles\tTemp(C)\tCORE_PERF_LIMIT_REASONS\tinstret\tllcmiss\tresource-stall\texec-stall\tbranch-retired\n");
		}
		unsigned long i;
		for (i = 0; i < runtime->sampler->samplectrs[j]; i++)
		{
			double time = (thread_samples[j][i + 1].tsc_data - thread_samples[j][i].tsc_data) /
				(FRQ_AS_GHZ(((thread_samples[j][i].frq_data & 0xFFFFul) >> 8)) * 1000000000.0);

			unsigned long diff = (thread_samples[j][i + 1].energy_data -
				thread_samples[j][i].energy_data);


			fprintf(outfile[j], 
				"%f\t%llx\t%llu\t%lf\t%llu\t%llx\t%llu\t%llu\t%llu\t%llu\t%llu\n",
				FRQ_AS_GHZ(((thread_samples[j][i].frq_data & 0xFFFFul) >> 8)),
				(unsigned long long) (thread_samples[j][i].frq_data & 0xFFFFul),
				//(unsigned long long) (thread_samples[j][i].tsc_data),
				(unsigned long long) (thread_samples[j][i + 1].tsc_data - thread_samples[j][i].tsc_data),
				diff * runtime->sys->rapl_energy_unit / time,
				/* (unsigned long long) ((thread_samples[j][i + 1].rapl_throttled & 0xFFFFFFFF) - (thread_samples[j][i].rapl_throttled & 0xFFFFFFFF)), */
				(unsigned long long) (80 - ((thread_samples[j][i].therm & 0x7F0000) >> 16)),
				(unsigned long long) thread_samples[j][i].perflimit,
				(unsigned long long) thread_samples[j][i + 1].instret - thread_samples[j][i].instret,
				(unsigned long long) thread_samples[j][i + 1].llcmiss - thread_samples[j][i].llcmiss,
				(unsigned long long) thread_samples[j][i + 1].restalls - thread_samples[j][i].restalls,
				(unsigned long long) thread_samples[j][i + 1].exstalls - thread_samples[j][i].exstalls,
				(unsigned long long) thread_samples[j][i + 1].branchret - thread_samples[j][i].branchret
				);
		}
	}
	last_dump_loc = runtime->sampler->samplectrs[0];
}

// display the usage of this program
void dump_help()
{
	//char validargs[] = "rlwLtdsveRhSAMO";
	printf("Valid options:\n");
	printf("\t-r: polling rate in samples per second\n");
	printf("\t-l: RAPL limit\n");
	printf("\t-L: hard RAPL limit\n");
	printf("\t-w: RAPL time window\n");
	printf("\t-t: clustering sensitivity threshold\n");
	printf("\t-d: display cutoff for profiles, if argument is missing then profiles are not dumped\n");
	printf("\t-s: system control, do sampling and profiling only\n");
	printf("\t-v: verbose\n");
	printf("\t-e: experimental feature enable\n");
	printf("\t-R: enable reporting\n");
	printf("\t-h: display this menu\n");
	printf("\t-S: disable memory power shifting\n");
	printf("\t-A: disable throttle avoidance\n");
	printf("\t-O: disable overpower\n");
	printf("\t-M: manual frequency override for memory phase (requires -S)\n");
	printf("\t-C: manual frequency override for compute phase (requires -A)\n");
}

// dump the configuration settings of this instance
void dump_config(struct powgov_runtime *runtime, FILE *out)
{
	fprintf(out, "Power Governor Configuration:\n");
	fprintf(out, "\trapl limit 1 %lf\n\trapl limit 2 %lf\n\trapl limit 1 time window %lf\n", 
			runtime->power->rapl1, runtime->power->rapl2, runtime->power->window);
	fprintf(out, "\tsamples per second %u\n", runtime->sampler->sps);
	fprintf(out, "\tgovernor bound to core %d\n", runtime->sys->num_cpu);
	if (runtime->cfg->man_cpu_ctrl == 1)
	{
		fprintf(out, "\tgovernor frequency control ON\n");
	}
	else
	{
		fprintf(out, "\tgovernor frequency control OFF\n");
	}
	if (runtime->cfg->mem_frq_override != 0.0)
	{
		fprintf(out, "\tforced memory phase frequency %lf\n", runtime->cfg->mem_frq_override);
	}
	if (runtime->cfg->cpu_frq_override != 0.0)
	{
		fprintf(out, "\tforced cpu phase frequency %lf\n", runtime->cfg->cpu_frq_override);
	}
	fprintf(out, "\tthrottle avoidance %s\n", (runtime->cfg->throttle_avoid ? "enabled" : "disabled"));
	fprintf(out, "\tmemory power shifting %s\n", (runtime->cfg->mem_pow_shift ? "enabled" : "disabled"));
	fprintf(out, "\tclassification threshold %lf\n", runtime->classifier->dist_thresh);
}

// dump system info
void dump_sys(struct powgov_runtime *runtime, FILE *out)
{
	fprintf(out, "System Configuration:\n");
	fprintf(out, "\tMax frequency %f\n\tBase Frequency %f\n\tSockets %d\n\tCores Per Socket %d\n",
			FRQ_AS_GHZ(runtime->sys->max_pstate), FRQ_AS_GHZ(runtime->sys->max_non_turbo),
			runtime->sys->sockets, runtime->sys->coresPerSocket);
	fprintf(out, "\tHyperthreads %s\n\tTotal Processors %d\n\tTDP %lf\n", 
			(runtime->sys->threadsPerCore ? "enabled" : "disabled"), runtime->sys->num_cpu, 
			runtime->power->proc_tdp);
	dump_rapl(out);
	hwpstuff(out);
	uint64_t misc = 0;
	read_msr_by_coord(0, 0, 0, 0x1A0, &misc);
	fprintf(out, "\tmisc enable %lx\n", misc);
	
}

// signal handler for power-governor, will shutdown and dump when it receives a SIGUSR1
void signal_exit(int signum)
{
	fprintf(stderr, "Sampler terminating...\n");
	fflush(stderr);
	LOOP_CTRL = 0;
	return;
}


int main(int argc, char **argv)
{
	// using libmsr
	if (init_msr())
	{
		fprintf(stderr, "ERROR: unable to init libmsr\n");
		exit(-1);
	}


	// have this process listen for SIGUSR1 signal
	struct sigaction sighand;
	memset(&sighand, 0, sizeof(struct sigaction));
	sighand.sa_handler = signal_exit;
	sigaction(SIGUSR1, &sighand, NULL);
	sigset_t sset;
	sigemptyset(&sset);
	sigaddset(&sset, SIGUSR1);
	sigprocmask(SIG_UNBLOCK, &sset, NULL);

	// initialize runtime configuration
	struct powgov_runtime *runtime = (struct powgov_runtime *) calloc(1,
			sizeof(struct powgov_runtime));
	runtime->cfg = (struct powgov_config *) calloc(1, sizeof(struct powgov_config));
	runtime->cfg->man_cpu_ctrl = 1; // -s for "system control"
	runtime->cfg->threadcount = 1;
	runtime->cfg->experimental = 0; // -e for experimental
	runtime->cfg->mem_pow_shift = 1;
	runtime->cfg->fup_timeout = 100;
	runtime->cfg->frq_duty_length = 10;
	runtime->cfg->cpu_frq_override = 0.0;
	runtime->cfg->mem_frq_override = 0.0;
	runtime->cfg->throttle_step_mult = 5.0;
	runtime->cfg->frq_change_step = 0.1;
	runtime->cfg->throttle_avoid = 1;

	// allocate all of the runtime data
	runtime->sys = (struct powgov_sysconfig *) calloc(1, sizeof(struct powgov_sysconfig));
	runtime->files = (struct powgov_files *) calloc(1, sizeof(struct powgov_files));
	runtime->sampler = (struct powgov_sampler *) calloc(1, sizeof(struct powgov_sampler));
	runtime->sampler->samplectrs = (unsigned long *) calloc(runtime->cfg->threadcount,
			sizeof(unsigned long));
	runtime->sampler->l1 = (struct powgov_l1 *) calloc(1, sizeof(struct powgov_l1));
	runtime->sampler->l2 = (struct powgov_l2 *) calloc(1, sizeof(struct powgov_l2));
	runtime->sampler->l3 = (struct powgov_l3 *) calloc(1, sizeof(struct powgov_l3));
	memset(runtime->sampler->l3->graph, 0, MAX_L3_GRAPH * sizeof(struct l3_graph_node));
	runtime->classifier = (struct powgov_classifier *) calloc(1,
			sizeof(struct powgov_classifier));
	runtime->power = (struct powgov_power *) calloc(1, sizeof(struct powgov_power));
	// initialize power governor configurization
	runtime->sampler->total_samples = 0;
	runtime->sampler->sps = 500; // -r for "rate"
	runtime->classifier->dist_thresh = 0.25;
	runtime->classifier->pct_thresh = 0.01;
	runtime->classifier->phase_cycle_thresh = 25000000;
	
	// TODO: sampling should just dump every x seconds
	VERBOSE = 0; // -v for "verbose"
	REPORT = 0; // -R for report

	// TODO: these should be read in from a file
	// these are the values at 800MHz, linear regression model based on frequency is used
	// cpu phase
	runtime->classifier->prof_class[0].ipc = 0.576;
	runtime->classifier->prof_class[0].mpc = 0.005;
	runtime->classifier->prof_class[0].rpc = 0.017;
	runtime->classifier->prof_class[0].epc = 0.027;
	runtime->classifier->prof_class[0].bpc = 0.0006;
	// memory phase
	runtime->classifier->prof_class[1].ipc = 0.122;
	runtime->classifier->prof_class[1].mpc = 0.004;
	runtime->classifier->prof_class[1].rpc = 0.118;
	runtime->classifier->prof_class[1].epc = 0.427;
	runtime->classifier->prof_class[1].bpc = 0.017;
	// IO/sleep phase (derived)
	runtime->classifier->prof_class[2].ipc = 0;
	runtime->classifier->prof_class[2].mpc = 0;
	runtime->classifier->prof_class[2].rpc = 0;
	runtime->classifier->prof_class[2].epc = 0;
	runtime->classifier->prof_class[2].bpc = 0;
	// mixed phase (derived)
	runtime->classifier->prof_class[3].ipc = 0.349;
	runtime->classifier->prof_class[3].mpc = 0.005;
	runtime->classifier->prof_class[3].rpc = 0.06;
	runtime->classifier->prof_class[3].epc = 0.25;
	runtime->classifier->prof_class[3].bpc = 0.005;
	// TODO: temporary solution
	frequency_scale_phase(&runtime->classifier->prof_class[0], 0.8, 4.5, 
			&runtime->classifier->prof_maximums);


	// lookup processor information with CPUID
	uint64_t rax, rbx, rcx, rdx;
	rax = rbx = rcx = rdx = 0;
	cpuid(0x16, &rax, &rbx, &rcx, &rdx);
	runtime->sys->min_pstate = 8;
	runtime->sys->max_pstate = ((rbx & 0xFFFFul) / 100);
	runtime->sys->max_non_turbo = ((rax & 0xFFFFul) / 100);
	runtime->classifier->mem_freq_throttle = (double) runtime->sys->max_pstate;
	uint64_t coresPerSocket, hyperThreads, sockets;
	int HTenabled;
	coresPerSocket = hyperThreads = sockets = HTenabled = 0;
	cpuid_detect_core_conf(&coresPerSocket, &hyperThreads, &sockets, &HTenabled);
	// TODO: currently hyper threads are ignored
	runtime->sys->num_cpu = sockets * coresPerSocket;
	assert(runtime->sys->num_cpu > 0);
	assert(runtime->sys->max_pstate > runtime->sys->max_non_turbo);
	assert(runtime->sys->max_non_turbo > runtime->sys->min_pstate);
	runtime->sys->sockets = sockets;
	runtime->sys->coresPerSocket = coresPerSocket;
	runtime->sys->threadsPerCore = hyperThreads;


	// process command line arguments
	// TODO this could be cleaner
	char validargs[] = "rlwLtdsveRhSAMOC";
	unsigned char numargs = strlen(validargs);
	int j;
	for (j = 1; j < argc; j++)
	{
		if (argv[j][0] != '-')
		{
			continue;
		}
		int badarg = 1;
		int i;
		for (i = 0; i < numargs; i++)
		{
			if (argv[j][1] == validargs[i])
			{
				badarg = 0;
				switch (argv[j][1])
				{
					case 'r':
						ARG_ERROR;
						runtime->sampler->sps = (unsigned) atoi(argv[j+1]);
						break;
					case 'l':
						ARG_ERROR;
						runtime->power->rapl1 = (double) atof(argv[j+1]);
						break;
					case 'w':
						ARG_ERROR;
						runtime->power->window = (double) atof(argv[j+1]);
						break;
					case 'L':
						ARG_ERROR;
						runtime->power->rapl2 = (double) atof(argv[j+1]);
						break;
					case 't':
						ARG_ERROR;
						runtime->classifier->dist_thresh = (double) atof(argv[j+1]);
						break;
					case 'd':
						ARG_ERROR;
						runtime->classifier->pct_thresh = (double) atof(argv[j+1]);
						break;
					case 's':
						runtime->cfg->man_cpu_ctrl = 0;	
						break;
					case 'v':
						VERBOSE = 1;
						break;
					case 'R':
						REPORT = 1;
						break;
					case 'e':
						runtime->cfg->experimental = 1;
						break;
					case 'O':
						printf("Error: this feature %s not implemented yet\n", argv[j]);
						exit(-1);
						break;
					case 'S':
						runtime->cfg->mem_pow_shift = 0;
						break;
					case 'A':
						runtime->cfg->throttle_avoid = 0;
						break;
					case 'M':
						ARG_ERROR;
						if (runtime->cfg->mem_pow_shift != 0)
						{
							printf("Error: %s requires memory power shifting to be disabled (-S)\n",
									argv[j]);
							exit(-1);
						}
						runtime->cfg->mem_frq_override = (double) atof(argv[j+1]) * 10;
						break;
					case 'C':
						ARG_ERROR;
						if (runtime->cfg->throttle_avoid != 0)
						{
							printf("Error: %s requires throttle avoidance to be disabled (-A)\n",
									argv[j]);
							exit(-1);
						}
						runtime->cfg->cpu_frq_override = (double) atof(argv[j+1]) * 10;
						break;
					case 'h':
					default:
						dump_help();
						exit(-1);
						break;
				}
			}
		}
		if (badarg)
		{
			printf("Error: invalid option %s\n", argv[j]);
			exit(-1);
		}
	}


	// finish configuration based on arguments
	unsigned srate = (1000.0 / runtime->sampler->sps) * 1000u;
	if (runtime->cfg->man_cpu_ctrl)
	{
		enable_turbo(runtime);
		set_perf(runtime, runtime->sys->max_pstate);
	}
	activate_performance_counters(runtime);


	// bind the power governor to core X
	cpu_set_t cpus;
	CPU_ZERO(&cpus);
	// TODO: make this an option
	CPU_SET(0, &cpus);
	sched_setaffinity(0, sizeof(cpus), &cpus);


	// setup RAPL
	uint64_t unit;
	read_msr_by_coord(0, 0, 0, MSR_RAPL_POWER_UNIT, &unit);
	uint64_t power_unit = unit & 0xF;
	runtime->sys->rapl_power_unit = 1.0 / (0x1 << power_unit);
	uint64_t seconds_unit_raw = (unit >> 16) & 0x1F;
	runtime->sys->rapl_seconds_unit = 1.0 / (0x1 << seconds_unit_raw);
	// TODO: loss of precision with energy unit
	uint64_t energy_unit_raw = ((unit >> 8) & 0x1F);
	runtime->sys->rapl_energy_unit = 1.0 / (0x1 << energy_unit_raw);
	uint64_t powinfo = 0;
	read_msr_by_coord(0, 0, 0, MSR_PKG_POWER_INFO, &powinfo);
	runtime->power->proc_tdp = (powinfo & 0x3FFF) * runtime->sys->rapl_power_unit;
	if (runtime->power->rapl1 == 0.0)
	{
		runtime->power->rapl1 = runtime->power->proc_tdp;
	}
	if (runtime->power->rapl2 == 0.0)
	{
		runtime->power->rapl2 = runtime->power->proc_tdp * 1.2;
	}
	set_rapl(runtime->power->window, runtime->power->rapl1, runtime->sys->rapl_power_unit, 
			runtime->sys->rapl_seconds_unit, 0);
	set_rapl2(100, runtime->power->rapl2, runtime->sys->rapl_power_unit, runtime->sys->rapl_seconds_unit, 0);
	runtime->power->energy_overflow = 0;


	// setup sampler levels
	runtime->sampler->l1->interval = 1000 / runtime->sampler->sps;
	runtime->sampler->l2->interval = (unsigned) (runtime->power->window * 1000.0);
	runtime->sampler->l3->interval = runtime->sampler->sps;
	runtime->sampler->l3->baseline_ipc = 0.0;
	runtime->power->excursion = 0;


	// print verbose descriptions to stdout
	if (VERBOSE)
	{
		printf(HEADERSTRING "\n\tPower Governor v%s\n" HEADERSTRING "\n", VERSIONSTRING);
		dump_sys(runtime, stdout);
		dump_config(runtime, stdout);
	}


	// begin sampling if argument present
	if (runtime->cfg->experimental)
	{
		init_sampling(runtime);
	}


	if (VERBOSE)
	{
		fprintf(stdout, "Initialization complete...\n");
	}


	// gather initial measurements for report if argument is present
	uint64_t inst_before, inst_after;
	uint64_t aperf_before, mperf_before;
	uint64_t ovf_stat;
	unsigned ovf_ctr = 0;
	double avgrate = 0.0;
	struct timeval start, current;
	uint64_t busy_tsc_pre[runtime->sys->num_cpu * 2]; // this only works because c99
	if (REPORT)
	{
		read_msr_by_coord(0, 0, 0, IA32_FIXED_CTR0, &inst_before);
		read_msr_by_coord(0, 0, 0, MSR_IA32_APERF, &aperf_before);
		read_msr_by_coord(0, 0, 0, MSR_IA32_MPERF, &mperf_before);
		//uint64_t busy_unh_pre[runtime->sys->num_cpu * 2];
		int j;
		for (j = 0; j < runtime->sys->num_cpu; j++)
		{
			read_msr_by_coord(0, j, 0, IA32_TIME_STAMP_COUNTER, &busy_tsc_pre[j]);
			read_msr_by_coord(1, j, 0, IA32_TIME_STAMP_COUNTER, &busy_tsc_pre[j + runtime->sys->num_cpu]);
			//read_msr_by_coord(0, j, 0, IA32_FIXED_CTR1, &busy_unh_pre[j]);
			//read_msr_by_coord(1, j, 0, IA32_FIXED_CTR1, &busy_unh_pre[j + runtime->sys->num_cpu]);
			write_msr_by_coord(0, j , 0, IA32_FIXED_CTR1, 0x0ul);
			write_msr_by_coord(1, j , 0, IA32_FIXED_CTR1, 0x0ul);
		}
	
	}


	// the main loop, continues until a SIGUSR1 is recieved
	uint64_t ovf_ctrl;
	gettimeofday(&start, NULL);
	sample_data(runtime);
	usleep(srate);
	while (LOOP_CTRL)
	{
		sample_data(runtime);
		read_msr_by_coord(0, 0, 0, IA32_PERF_GLOBAL_STATUS, &ovf_stat);
		if (ovf_stat & 0x1)
		{
			read_msr_by_coord(0, 0, 0, IA32_PERF_GLOBAL_OVF_CTRL, &ovf_ctrl);
			write_msr_by_coord(0, 0, 0, IA32_PERF_GLOBAL_OVF_CTRL, ovf_ctrl & 0xFFFFFFFFFFFFFFFE);
			ovf_ctr++;
		}
		l1_analysis(runtime);
		usleep(srate);
	}
	gettimeofday(&current, NULL);
	double exec_time = (double) (current.tv_sec - start.tv_sec) +
			(current.tv_usec - start.tv_usec) / 1000000.0;


	if (VERBOSE)
	{
		fprintf(stdout, "Power Governor terminating...\n");
	}


	// dump all report files if argument present
	if (REPORT)
	{
		uint64_t aperf_after, mperf_after;
		read_msr_by_coord(0, 0, 0, MSR_IA32_APERF, &aperf_after);
		read_msr_by_coord(0, 0, 0, MSR_IA32_MPERF, &mperf_after);
		read_msr_by_coord(0, 0, 0, IA32_FIXED_CTR0, &inst_after);

		avgrate = exec_time / runtime->sampler->total_samples;
		runtime->files->sreport = fopen("sreport", "w");
		dump_sys(runtime, runtime->files->sreport);
		dump_config(runtime, runtime->files->sreport);
		fprintf(runtime->files->sreport, "Results:\n");
		fprintf(runtime->files->sreport, "\tActual run time: %f\n", exec_time);
		fprintf(runtime->files->sreport, "\tAverage Sampling Rate: %lf seconds\n", avgrate);
		fprintf(runtime->files->sreport, "\tAvg Frq: %f\n", (float) 
				(aperf_after - aperf_before) / (float) 
				(mperf_after - mperf_before) * FRQ_AS_GHZ(runtime->sys->max_non_turbo));
		fprintf(runtime->files->sreport, "\tInstructions: %lu (ovf %u)\n", 
				inst_after - inst_before, ovf_ctr);
		fprintf(runtime->files->sreport, "\tIPS: %lf (ovf %u)\n", 
				(inst_after - inst_before) / exec_time, ovf_ctr);
		
		uint64_t energy_end = runtime->sampler->l1->new_sample.energy_data;
		uint64_t energy_begin = runtime->sampler->first_sample.energy_data;
		if (runtime->power->energy_overflow > 0)
		{
			fprintf(runtime->files->sreport, "\tTotal Power: %lf\n", (double) ((0xEFFFFFFF * runtime->power->energy_overflow - 
					energy_begin + energy_end) * runtime->sys->rapl_energy_unit) / exec_time);
			fprintf(runtime->files->sreport, "\tEnergy Overflows: %u\n", runtime->power->energy_overflow);
		}
		else
		{
			fprintf(runtime->files->sreport, "\tTotal Power: %lf\n", ((energy_end - energy_begin) * 
						runtime->sys->rapl_energy_unit) / exec_time);
		}

		char fname[FNAMESIZE];
		snprintf((char *) fname, FNAMESIZE, "powgov_profiles");
		runtime->files->profout = fopen(fname, "w");
		// TODO: figure out if miscounts from remove_unused or something else (fixed?)
		//remove_unused(runtime);
		//agglomerate_profiles(runtime);
		dump_phaseinfo(runtime, runtime->files->profout, &avgrate);
		int j;
		for (j = 0; j < runtime->sys->num_cpu; j++)
		{
			uint64_t busy_tsc_post, busy_unh_post;
			read_msr_by_coord(0, j, 0, IA32_TIME_STAMP_COUNTER, &busy_tsc_post);
			read_msr_by_coord(0, j, 0, IA32_FIXED_CTR1, &busy_unh_post);
			double diff_unh = ((double) busy_unh_post);
			double diff_tsc = ((double) busy_tsc_post - busy_tsc_pre[j]);
			fprintf(runtime->files->sreport, "\ts1c%d pct busy: %lf%% (%lf/%lf)\n", j,  
					diff_unh / diff_tsc * 100.0, diff_unh, diff_tsc);
			//read_msr_by_coord(1, j, 0, IA32_TIME_STAMP_COUNTER, &busy_tsc_post);
			//read_msr_by_coord(1, j, 0, IA32_FIXED_CTR1, &busy_unh_post);
			//diff_unh = ((double) busy_unh_post);
			//diff_tsc = ((double) busy_tsc_post - busy_tsc_pre[j + runtime->sys->num_cpu]);
			//printf("s2c%d pct busy: %lf\% (%lf/%lf)\n", j,  diff_unh / diff_tsc * 100.0, diff_unh, diff_tsc);
		}
		fclose(runtime->files->sreport);
		fclose(runtime->files->profout);
	}


	// dump sampler data if argument present
	if (runtime->cfg->experimental)
	{
		dump_data(runtime, runtime->files->sampler_dumpfiles);
		int i;
		for (i = 0; i < runtime->cfg->threadcount; i++)
		{
			free(runtime->sampler->thread_samples[i]);
			fclose(runtime->files->sampler_dumpfiles[i]);
		}
		free(runtime->sampler->thread_samples);
		free(runtime->sampler->samplectrs);
		free(runtime->files->sampler_dumpfiles);
	}


	// Reset RAPL to defaults
	set_rapl(1, runtime->power->proc_tdp, runtime->sys->rapl_power_unit, runtime->sys->rapl_seconds_unit, 0);
	set_rapl2(100, runtime->power->proc_tdp * 1.2, runtime->sys->rapl_power_unit, 
			runtime->sys->rapl_seconds_unit, 0);
	finalize_msr();
	return 0;
}
