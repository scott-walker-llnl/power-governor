#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "powgov_l3.h"
#include "powgov_l1.h"
#include "powgov.h"

// "freeze" a phase, which prevents the level 3 CPU frequency from being changed by level 1
static void l3_phase_freeze(struct powgov_runtime *runtime, struct l3_graph_node *node)
{
	int i;
	// find the minimum throttle freq
	for (
		i = 0;
		i < MAX_PSTATE_NUM && node->next[0]->phase->workload.throttlehist[i] == 0;
		i++);
	if (i == MAX_PSTATE_NUM)
	{
		// this cpu phase has not been throttled, there is no benefit in
		// throttling the memory phase
		return;
	}

	printf("\tmin throttle freq is %f (count %d)\n", FRQ_AS_GHZ(i + 8),
			node->next[0]->phase->workload.throttlehist[i]);

	// freeze at this frequency
	node->next[0]->phase->l3_freeze = i;

	printf("\tcpu phase %d frozen at %f\n", (int)(node->next[0]->phase -
			runtime->classifier->phases), FRQ_AS_GHZ(i + 8));

	node->phase->workload.frq_target = runtime->sys->max_pstate - 7.0;
	node->phase->workload.flags |= L3_FREEZE_FLAG;
	node->next[0]->phase->workload.flags |= L3_FREEZE_FLAG;

	memcpy(node->phase->pthrottlehist, node->phase->workload.throttlehist, sizeof(short) *
			MAX_PSTATE_NUM);
	memcpy(node->next[0]->phase->pthrottlehist, node->next[0]->phase->workload.throttlehist,
			sizeof(short) * MAX_PSTATE_NUM);
	memset(node->phase->workload.throttlehist, 0, sizeof(short) * MAX_PSTATE_NUM);
	memset(node->next[0]->phase->workload.throttlehist, 0, sizeof(short) * MAX_PSTATE_NUM);

	printf("\tphase %d freq set to %f\n", (int) (node->phase - runtime->classifier->phases), FRQ_AS_GHZ(node->phase->workload.frq_target));
	// TODO: should we make the current workload use this frequency if it is this phase?
}

// check a frozen phase to see if something good or bad happened
static void l3_freeze_report(struct powgov_runtime *runtime, struct l3_graph_node *node)
{
	int minfrqidx = node->next[0]->phase->l3_freeze;
	printf("\tthrottles at frq %f: before %d after %d\n", FRQ_AS_GHZ(minfrqidx + 8),
			node->next[0]->phase->pthrottlehist[minfrqidx],
			node->next[0]->phase->workload.throttlehist[minfrqidx]);

	if (node->next[0]->phase->pthrottlehist[minfrqidx] >
			node->next[0]->phase->workload.throttlehist[minfrqidx])
	{
		// throttling helped
		printf("\tthrottling phase %d was good\n", (int)(node->phase - runtime->classifier->phases));
	}
	else
	{
		// throttling didn't help
		printf("\tthrottling phase %d was bad\n", (int)(node->phase - runtime->classifier->phases));
		// unfreeze mem phase
		node->phase->workload.flags &= ~L3_FREEZE_FLAG;
		node->phase->workload.frq_target = runtime->sys->max_pstate;
	}
	memset(node->phase->workload.throttlehist, 0, sizeof(short) * MAX_PSTATE_NUM);
	memset(node->next[0]->phase->workload.throttlehist, 0, sizeof(short) * MAX_PSTATE_NUM);
	// unfreeze the cpu phase
	node->next[0]->phase->workload.flags &= ~L3_FREEZE_FLAG;
}

// look for MEM phases occurring before CPU phases, attempt MEM throttling for power savings
void l3_analysis(struct powgov_runtime *runtime)
{
	// look through graph and identify MEM phases before CPU phases
	// TODO: check if end-to-end power over this window is below limit
	printf("\nl3 analysis\n");
	if (runtime->cfg->mem_pow_shift == 0)
	{
		printf("mem pow shift disabled\n");
		return;
	}
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH; itr++)
	{
		if (itr->phase == NULL || itr->next[0] == NULL)
		{
			continue;
		}
		// TODO: if phase didn't occur, we don't want to do anything at all
		if (itr->phase->workload.class == CLASS_MEM)
		{
			if (itr->next[0]->phase->workload.class == CLASS_CPU)
			{
				printf("\tphase %d is a mem before cpu phase\n", (int) (itr->phase - runtime->classifier->phases));
				// we have found a mem phase that occurs before a cpu phase,
				// we should freeze the cpu phase to its minimum throttle freq
				// then next time it comes around we can check to see if we reduced
				// the number of throttles that occurs at minimum throttle freq
				if (itr->phase->l3_freeze == 0) 
				{
					if ((itr->phase->workload.flags & L3_FREEZE_FLAG) &&
						(itr->next[0]->phase->workload.flags & L3_FREEZE_FLAG))
					{
						l3_freeze_report(runtime, itr);
					}
					else
					{
						l3_phase_freeze(runtime, itr);
					}
					// we always want to freeze this phase to postpone re-analysis
					itr->phase->l3_freeze = L3_FREEZE_LEN;
				}
				itr->phase->l3_freeze--;
			}
		}
	}
}

// decay the count of successor occurrences to kick out stale entries
static inline void decay_occurrence(struct powgov_runtime *runtime)
{
	struct l3_graph_node *current = runtime->sampler->l3->current_node;
	int i;
	for (i = 0; i < 3; i++)
	{
		if (current->next_occurrences[i] >= 1)
		{
			current->next_occurrences[i]--;
		}
		if (current->next_occurrences[i] == 0 && current->next[i] != NULL)
		{
			// delete this entry and slide things over
			if (i < 2)
			{
				if (current->next[i]->reference_counter >= 1)
				{
					current->next[i]->reference_counter--;
				}
				current->next[i] = current->next[i + 1];
				current->next_occurrences[i] = current->next_occurrences[i + 1];
				current->next[i + 1] = NULL;
				current->next_occurrences[i + 1] = 0;
				// skip the next iteration since that one is gone
				i++;
			}
		}
	}
}

// update the last entry (or first NULL) in the successor list
static inline void update_last(struct powgov_runtime *runtime, struct l3_graph_node *node)
{
	struct l3_graph_node *current = runtime->sampler->l3->current_node;
	int i;
	// get either unused position or last one
	for (i = 0; i < 2 && current->next[i] != NULL; i++);
	current->next[i] = node;
	// seed this with a large enough number that it doesn't get evicted immediately
	current->next_occurrences[i] = 30;
}

// do a simple sort so that the first successor is always the most occurring one
static inline void sort_by_occurrence(struct powgov_runtime *runtime)
{
	struct l3_graph_node *current = runtime->sampler->l3->current_node;
	int i;
	/* printf("before sort node %d\n", (int) (current->phase - runtime->classifier->phases)); */
	/* for (i = 0; i < 3; i++) */
	/* { */
	/* 	if (current->next[i] != NULL) */
	/* 	{ */
	/* 		printf("\tsuccessor %d count %lu\n", i, current->next_occurrences[i]); */
	/* 	} */
	/* } */
	for (i = 2; i > 0; i--)
	{
		// swap left while occurrences is greater
		if (current->next[i] != NULL &&
			current->next_occurrences[i - 1] < current->next_occurrences[i])
		{
			uint64_t temp_occ = current->next_occurrences[i - 1];
			struct l3_graph_node *temp_next = current->next[i - 1];
			current->next_occurrences[i - 1] = current->next_occurrences[i];
			current->next[i - 1] = current->next[i];
			current->next_occurrences[i] = temp_occ;
			current->next[i] = temp_next;
		}
	}
	/* printf("after sort node %d\n", (int) (current->phase - runtime->classifier->phases)); */
	/* for (i = 0; i < 3; i++) */
	/* { */
	/* 	if (current->next[i] != NULL) */
	/* 	{ */
	/* 		printf("\tsuccessor %d count %lu\n", i, current->next_occurrences[i]); */
	/* 	} */
	/* } */
}

// add a new graph node when we find a new phase
void add_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase)
{
	// initialize a new graph node object pointing to new_phase
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != NULL; itr++);
	itr->phase = new_phase;
	itr->reference_counter = 0;
	memset(itr->next, 0, 3 * sizeof(struct l3_graph_node *));
	memset(itr->next_occurrences, 0, 3 * sizeof(uint64_t));
	// set current->next to the new graph node
	if (runtime->sampler->l3->current_node != NULL)
	{
		itr->reference_counter++;
		update_last(runtime, itr);
		sort_by_occurrence(runtime);
		decay_occurrence(runtime);
	}
	// set current = newly added node
	runtime->sampler->l3->current_node = itr;
}

// if there is already a successor matching this phase, update its successor count
static void update_existing_next(struct powgov_runtime *runtime, struct l3_graph_node *node)
{
	struct l3_graph_node *current = runtime->sampler->l3->current_node;
	int i;
	// find existing or place it at end
	for (i = 0; i < 2 && current->next[i] != node; i++);
	if (i == 2)
	{
		// we are replacing the end
		current->next[i] = node;
		node->reference_counter++;
		current->next_occurrences[i] = 2;
	}
	else
	{
		// we are updating an existing
		current->next_occurrences[i] += 2;
	}
	sort_by_occurrence(runtime);
	decay_occurrence(runtime);
}

// find an existing graph node and update its successor list to point to the found phase
void update_graph_node(struct powgov_runtime *runtime, struct phase_profile *updated_phase)
{
	// find the existing graph node
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = runtime->sampler->l3->graph;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != updated_phase; itr++);

	// set current->next = found node
	update_existing_next(runtime, itr);

	// set current = found node
	runtime->sampler->l3->current_node = itr;
}
