#include <string.h>
#include "powgov_l3.h"
#include "powgov_l1.h"

void l3_analysis(struct powgov_runtime *runtime)
{
	// check how ipc of various phases changes over time in response to L1
}

struct phase_profile *update_phase(struct powgov_runtime *runtime, struct workload_profile *this_profile, struct phase_profile *prof, double phase_cycles)
{
	prof->phase_occurrences++;
	uint64_t prev_occurrences = prof->workload.occurrences;
	uint64_t new_occurrences = this_profile->occurrences;

	prof->workload.occurrences = prev_occurrences + new_occurrences;

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

void add_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase)
{
	// initialize a new graph node object pointing to new_phase
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != NULL; itr++);
	// initialize new graph node -> nextnode = NULL
	itr->nextnode = NULL;
	itr->phase = new_phase;
	// set predicted = NULL because we haven't seen this before
	// TODO: make sure predicted phase is not used iff NULL
	runtime->sampler->l3->predicted_phase = NULL;
	// set current->next to the new graph node
	if (runtime->sampler->l3->current_node != NULL)
	{
		runtime->sampler->l3->current_node->nextnode = itr;
	}
	// set current = newly added node
	runtime->sampler->l3->current_node = itr;
}

void update_graph_node(struct powgov_runtime *runtime, struct phase_profile *updated_phase)
{
	// find the existing graph node
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != updated_phase; itr++);
	// set current->next = found node
	runtime->sampler->l3->current_node->nextnode = itr;
	// set prediction = found node -> next -> phase
	runtime->sampler->l3->predicted_phase = itr->nextnode->phase;
	// set current = found node
	runtime->sampler->l3->current_node = itr;
}
