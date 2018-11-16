#pragma once
#include <stdint.h>
#include "powgov.h"
#include "powgov_profiles.h"

#define MAX_L3_SEQ 32
#define MAX_L3_GRAPH 64
#define L3_FREEZE_LEN 4
#define L3_FREEZE_FLAG 0x1

struct powgov_runtime;

struct l3_graph_node
{
	struct phase_profile *phase; // the phase associated with this graph node
	struct l3_graph_node *next[3]; // the phases which succeed the current phase
	uint64_t next_occurrences[3]; // how many times the successors have occurred
	unsigned short reference_counter;
};

struct powgov_l3
{
	struct data_sample new_sample; // the new sample at the l3 sample rate
	struct data_sample prev_sample; // the previous sample at the l3 sample rate
	unsigned long interval; // the sample rate for this level
	uint64_t baseline_ipc; 
	//struct phase_profile sequence[MAX_L3_SEQ];
	struct l3_graph_node graph[MAX_L3_GRAPH]; // stores the graph nodes
	struct l3_graph_node *current_node; // the current location in the graph
	struct l3_graph_node *predictor; // the predicted location for the next phase
};

void l3_analysis(struct powgov_runtime *runtime);
void add_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase);
void update_graph_node(struct powgov_runtime *runtime, struct phase_profile *updated_phase);
