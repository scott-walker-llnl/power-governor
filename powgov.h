#pragma once
#include <stdio.h>
#include "msr_core.h"
#include "master.h"
#include "powgov_profiles.h"
#include "powgov_sampler.h"

#define FNAMESIZE 32

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
	double rapl1;
	double rapl_scaled;
	double rapl2;
	double window;
	double proc_tdp;
	uint64_t energy_begin;
	unsigned energy_overflow;
	char excursion;
};

struct powgov_config
{
	char man_cpu_ctrl;
	char throttle_avoid;
	char mem_pow_shift;
	char experimental;
	float cpu_frq_override;
	float mem_frq_override;
	float frq_change_step;
	float throttle_step_mult;
	unsigned short fup_timeout;
	unsigned short frq_duty_length;
	unsigned short threadcount;
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
