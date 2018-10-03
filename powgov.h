#pragma once
#include "msr_core.h"
#include "master.h"
#include "powgov_profiles.h"
#include "powgov_sampler.h"
#include "powgov_l1.h"
#include "powgov_l2.h"
#include "powgov_l3.h"

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

struct powgov_sampler
{
	unsigned long *samplectrs; // counter for number of samples on each thread
	struct data_sample **thread_samples;
	unsigned long numsamples;
	struct data_sample first_sample;
	struct powgov_l1 l1;
	struct powgov_l2 l2;
	struct powgov_l3 l3;
	unsigned sps;
	unsigned long total_samples;
};

struct powgov_runtime
{
	struct powgov_sysconfig *sys;
	struct powgov_files *files;
	struct powgov_sampler *sampler;
	struct powgov_classifier *classifier;
	struct powgov_power *power;
};

void dump_phaseinfo(struct powgov_runtime *runtime, FILE *outfile, double *avgrate);
void dump_data(struct powgov_runtime *runtime, FILE **outfile);
void dump_help();
void dump_config(struct powgov_runtime *runtime, FILE *out);
void dump_sys(struct powgov_runtime *runtime, FILE *out);
