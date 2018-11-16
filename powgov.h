#pragma once
#include <stdio.h>
#include "msr_core.h"
#include "master.h"
#include "powgov_profiles.h"
#include "powgov_sampler.h"

#define FNAMESIZE 32

#define FRQ_AS_GHZ(mfrq) ((mfrq) / 10.0)

struct powgov_classifier;
struct powgov_sampler;
struct data_sample;

struct powgov_sysconfig
{
	unsigned num_cpu; // number of cores present
	double rapl_energy_unit;
	double rapl_seconds_unit;
	double rapl_power_unit;
	double max_non_turbo;
	double max_pstate;
	double min_pstate;
	unsigned sockets;
	unsigned coresPerSocket;
	unsigned threadsPerCore;
};

struct powgov_files
{
	FILE *sreport; // main report file generated
	FILE **sampler_dumpfiles; // report files for sampler
	FILE *profout; // report file for profile clusters
};

struct powgov_power
{
	double rapl1; // the power limit for RAPL 1
	double rapl_scaled;
	double rapl2; // the power limit for RAPL 2
	double window; // the time window for RAPL 1
	double proc_tdp; // the TDP of the processor
	uint64_t energy_begin; // an energy sample taken at the beginning of the program
	unsigned energy_overflow; // count how many times the energy counter overflowed
	char excursion; // is true if the power limit has been exceeded
};

struct powgov_config
{
	char man_cpu_ctrl; // if true power-governor controlls processor frequency
	char throttle_avoid; // if true RAPL based throttles will be avoided
	char mem_pow_shift; // if true will attempt to save power in memory phases
	char experimental; // if true will enable experimental features
	float cpu_frq_override; // force CPU phases to a certain frequency
	float mem_frq_override; // force MEM phases to a certain frequency
	float frq_change_step; // how many MHz to change frequency by every iteration
	float throttle_step_mult; // multiplier for frq_change_step as multiplicative decrease
	unsigned short fup_timeout; // timeout in powgov samples before rasing phase freq
	unsigned short frq_duty_length; // used to duty cycle frequency for frq in 10s of MHz
	unsigned short threadcount; // how many threads to sample on (unsupported currently)
};

struct powgov_runtime
{
	struct powgov_sysconfig *sys;
	struct powgov_files *files;
	struct powgov_sampler *sampler;
	struct powgov_classifier *classifier;
	struct powgov_power *power;
	struct powgov_config *cfg;
};

void dump_phaseinfo(struct powgov_runtime *runtime, FILE *outfile, double *avgrate);
void dump_data(struct powgov_runtime *runtime, FILE **outfile);
void dump_help();
void dump_config(struct powgov_runtime *runtime, FILE *out);
void dump_sys(struct powgov_runtime *runtime, FILE *out);
void dump_error_report(struct powgov_runtime *runtime);
