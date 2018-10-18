#include <string.h>
#include <stdlib.h>
#include "powgov_l3.h"
#include "powgov_l1.h"

void l3_analysis(struct powgov_runtime *runtime)
{
	// look through graph and identify MEM phases before CPU phases
	// TODO: check if end-to-end power over this window is below limit
	printf("l3 analysis\n");
	if (runtime->cfg->mem_pow_shift == 0)
	{
		printf("mem pow shift disabled\n");
		return;
	}
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH; itr++)
	{
		if (itr->phase == NULL)
		{
			continue;
		}
		if (itr->phase->workload.class == CLASS_MEM && itr->nextnode != NULL)
		{
			if (itr->nextnode->phase->workload.class == CLASS_CPU)
			{
				if (itr->phase->l3_freeze == 0) 
				{
					printf("throttling mem phase\n");
					print_profile(&itr->phase->workload);
					// we have found a MEM before CPU phase, drop target freq
					itr->phase->workload.frq_target -= runtime->cfg->frq_change_step * 2.0;
					//  freeze this phase for a while so that we don't drop too fast
					itr->phase->l3_freeze = 3;
					itr->phase->ipc_history = itr->phase->workload.ipc;
					itr->nextnode->phase->ipc_history = itr->nextnode->phase->workload.ipc;
					itr->nextnode->phase->l3_freeze = 1;
				}
				else
				{
					// check if ipc increases for CPU phases after throttled MEM phases 
					double memdiff = itr->phase->workload.ipc - itr->phase->ipc_history;
					double cpudiff = itr->nextnode->phase->workload.ipc -
							itr->nextnode->phase->ipc_history;
					// if memdiff is negative, we slowed down memory phase
					// 		memdiff is >= 0 we did not effect memory phase
					// if cpudiff is negative, we slowed down cpu phase
					// 		cpudiff is >= 0 we sped up or did not effect cpu phase
					// if -mem, +cpu GOOD
					// if -mem, -cpu BAD
					// if +mem, +cpu GOOD
					// if +mem, -cpu BAD
					if (memdiff < 0.0 && cpudiff >= 0.0)
					{
						// this is good, mem was slowed but overall IPC increased
						// we at least saved power if broke even
						itr->phase->l3_freeze--;
					}
					else if (memdiff >= 0.0 && cpudiff >= 0.0)
					{
						// this is good, we either sped up one or both without hurting other
						itr->phase->l3_freeze--;
					}
					else
					{
						printf("unthrottling mem phase\n");
						print_profile(&itr->phase->workload);
						// either of following
						// (memdiff < 0.0 && cpudiff < 0.0)
						// (memdiff >= 0.0 && cpudif < 0.0)
						// this is bad, we made CPU or everything worse
						itr->phase->l3_freeze++;
						if (itr->phase->l3_freeze > 5 &&
							itr->phase->workload.frq_target < runtime->sys->max_pstate)
						{
							// this was consistently a bad decision and we are still under
							// frequency, so we should raise it
							itr->phase->l3_freeze = 3;
							itr->phase->workload.frq_target += 
									runtime->cfg->frq_change_step * 2.0;
						}
					}
				}
			}
		}
	}
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
