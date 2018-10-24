#pragma once
#include <stdint.h>
#include "powgov.h"
#include "powgov_profiles.h"

#define MAX_L3_SEQ 32
#define MAX_L3_GRAPH 128

struct powgov_runtime;

struct l3_graph_node
{
	struct phase_profile *phase;
	struct l3_graph_node *next; // most recently used successor
	char l3_freeze; // used to prevent l3 from dropping again too quickly
	char l3_thaw; // used to prevent l3 from raising too quickly
	char pow_shift; // if true this phase is involved in power shifting
	double ipc_history; // ipc before l3 intervention
};

struct powgov_l3
{
	struct data_sample new_sample;
	struct data_sample prev_sample;
	unsigned long interval;
	uint64_t baseline_ipc;
	//struct phase_profile sequence[MAX_L3_SEQ];
	struct l3_graph_node graph[MAX_L3_GRAPH];
	struct l3_graph_node *current_node;
	struct l3_graph_node *prelim;
};

void l3_analysis(struct powgov_runtime *runtime);
void add_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase);
void update_graph_node(struct powgov_runtime *runtime, struct phase_profile *updated_phase);
