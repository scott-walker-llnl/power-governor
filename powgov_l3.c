#include <string.h>
#include <stdlib.h>
#include "powgov_l3.h"
#include "powgov_l1.h"

// TODO: mem pow shift needs to have separate action for MRU and MOS
static void l3_mem_freeze(struct powgov_runtime *runtime, struct l3_graph_node *mem_node)
{
	printf("freezing mem phase\n");
	print_profile(&mem_node->phase->workload);

	if (mem_node->pow_shift == 0)
	{
		mem_node->ipc_history = mem_node->phase->workload.ipc;
		mem_node->next->ipc_history = mem_node->next->phase->workload.ipc;
	}

	mem_node->pow_shift = 1;
	mem_node->next->pow_shift = 1;

	if (mem_node->phase->workload.frq_target > runtime->sys->min_pstate)
	{
		mem_node->phase->workload.frq_target -= 1.0; // TODO: use config for this
	}
	// defer another frequency step until we've done analysis a few times 
	mem_node->l3_freeze = 3;
	mem_node->l3_thaw = 3;
}

static void l3_mem_thaw(struct powgov_runtime *runtime, struct l3_graph_node *mem_node)
{
	printf("thawing mem phase\n");
	print_profile(&mem_node->phase->workload);

	if (mem_node->phase->workload.frq_target < runtime->sys->max_pstate)
	{
		mem_node->phase->workload.frq_target += 1.0; // TODO: use config for this
	}
	else if (mem_node->phase->workload.frq_target == runtime->sys->max_pstate)
	{
		// we reached maximum frequency again, if we start shifting again then we need to
		// re-collect ipc history
		mem_node->pow_shift = 0;
		mem_node->next->pow_shift = 0;
	}

	// defer another frequency step until we've done analysis a few times 
	mem_node->l3_thaw = 3;
	mem_node->l3_freeze = 3;
}

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
		if (itr->phase->workload.class == CLASS_MEM && itr->next != NULL)
		{
			if (itr->next->phase->workload.class == CLASS_CPU)
			{
				// we have found a MEM node with a CPU successor in the graph
				if (itr->next->phase->workload.frq == runtime->sys->max_pstate)
				{
					// the CPU successor is already at maximum frequency so we can't help it
					continue;
				}
				if (itr->l3_freeze == 0.0)
				{
					l3_mem_freeze(runtime, itr);
					// we changed things so don't analyze this round
					return;
				}
				else if (itr->l3_thaw == 0.0)
				{
					l3_mem_thaw(runtime, itr);
					// we changed things so don't analyze this round
					return;
				}
				// check if ipc increases for CPU phases after throttled MEM phases 
				double memdiff = itr->phase->workload.ipc - itr->ipc_history;
				double cpudiff = itr->next->phase->workload.ipc -
						itr->next->ipc_history;
				// if -mem, +cpu GOOD
				// if -mem, -cpu BAD
				// if +mem, +cpu GOOD
				// if +mem, -cpu BAD
				if (cpudiff >= 0.0)
				{
					// we increased CPU ipc so maybe we helped
					double netipc = memdiff + cpudiff;
					if (netipc >= 0.0)
					{
						// we increased CPU ipc more than decreased MEM ipc
						// definitely helped
						itr->l3_freeze--;
						itr->l3_thaw++;
						// TODO: update history here?
					}
					else
					{
						// we increased CPU ipc less than decreased MEM ipc
						// throttled too much
						// don't delay more freezing, might be a one-off
						// encourage thawing because throttling likely too aggressive
						itr->l3_thaw--;
					}
				}
				else
				{
					// we didn't increase CPU ipc so didn't help
					// delay freezing because we don't want fequency to drop more
					// don't encourage thawing
					itr->l3_freeze++;
				}
			}
		}
	}
}

static struct l3_graph_node *find_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase)
{
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = begin;
	struct l3_graph_node *prelim = runtime->sampler->l3->prelim;
	// see if there exists a graph node that matches the current preliminary's phase and has
	// a next that matches the new phase
	for (; itr < begin + MAX_L3_GRAPH; itr++)
	{
		if (itr->phase == prelim->phase && itr->next != NULL && itr->next->phase == new_phase)
		{
			// a graph node for this parent/child combination already exists
			// point current to that instead of creating a new one
			return itr;
		}
	}
	return NULL;
}

static struct l3_graph_node *get_new_graph_node(struct powgov_runtime *runtime)
{
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = begin;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != NULL; itr++);
	if (itr == begin + MAX_L3_GRAPH && itr->phase != NULL)
	{
		printf("ERROR: out of graph storage\n");
		exit(-1);
	}
	return itr;
}

static void reject_preliminary(struct powgov_runtime *runtime, struct l3_graph_node *node)
{
	struct l3_graph_node *begin = runtime->sampler->l3->graph;
	struct l3_graph_node *itr = begin;
	for (; itr < begin + MAX_L3_GRAPH && itr->phase != NULL; itr++)
	{
		if (itr->next == runtime->sampler->l3->prelim)
		{
			// we don't wan't anything pointing to prelim because it was rejected
			// update everything pointing to prelim to point to the already seen node
			itr->next = node;
		}
	}

	// delete the current preliminary from the graph buffer
	memset(runtime->sampler->l3->prelim, 0, sizeof(struct l3_graph_node));

	// prelim was rejected, we have attached back to the graph
	runtime->sampler->l3->prelim = NULL;
}

static void dump_graph(struct powgov_runtime *runtime)
{
	/* if (runtime->sampler->l3->current_node->next_mos != itr) */
	/* { */
	/* 	printf("\nupdated phase ID %d to preceed phase ID %d\n",  */
	/* 			(int) (runtime->sampler->l3->current_node->phase - */
	/* 				runtime->classifier->phases), */
	/* 			(int) (updated_phase - runtime->classifier->phases)); */
	/* } */
	// set prediction = found node -> next -> phase
	/* printf("\tphase ID %d predicted next phase %d\n",  */
	/* 		(int) (itr->phase - runtime->classifier->phases), */
	/* 		(int) (runtime->sampler->l3->predicted_phase - runtime->classifier->phases)); */
	// set current = found node
	char loops[MAX_L3_GRAPH];
	char visited[MAX_L3_GRAPH];
	int i;
	for (i = 0; i < MAX_L3_GRAPH; i++)
	{
		if (runtime->sampler->l3->graph[i].phase == NULL)
		{
			visited[i] = 1;
		}
		else
		{
			visited[i] = 0;
		}
	}
	
	printf("\n");
	for (i = 0; i < MAX_L3_GRAPH; i++)
	{
		if (visited[i] == 0)
		{
			memset(loops, 0, MAX_L3_GRAPH);
			struct l3_graph_node *node = &runtime->sampler->l3->graph[i];
			printf("dumping graph\n");
			while (node != NULL && loops[node - runtime->sampler->l3->graph] == 0)
			{
				printf("%d -> ", (int) (node->phase - runtime->classifier->phases));
				visited[node - runtime->sampler->l3->graph] += 1;
				loops[node - runtime->sampler->l3->graph] += 1;
				node = node->next;
			}
			if (node != NULL && node->next != NULL && loops[node - runtime->sampler->l3->graph])
			{
				printf("(loops to %d)", node->next->phase - runtime->classifier->phases);
			}
			printf("\n");
		}
	}
}

// used to add node for phase never seen before
void add_graph_node(struct powgov_runtime *runtime, struct phase_profile *new_phase)
{
	// initialize a new graph node object pointing to new_phase
	struct l3_graph_node *new = get_new_graph_node(runtime);
	// initialize new graph node -> next... = NULL
	new->next = NULL;
	new->phase = new_phase;
	printf("new node\n");

	// always update current->next
	if (runtime->sampler->l3->current_node != NULL)
	{
		// if prelim exists, current should point to it, this will update its next
		// otherwise, current is a recently added node with no next
		// the preliminary must be added: since node is new preliminary it must be unique
		// we don't need to update pointers since the preliminary was not rejected
		runtime->sampler->l3->current_node->next = new;
	}

	// prelim either nonexistent or just committed
	runtime->sampler->l3->prelim = NULL;

	// set current = newly added node
	runtime->sampler->l3->current_node = new;

	dump_graph(runtime);
}

// we have seen this phase before, create preliminary because we don't know who comes next yet
void update_graph_node(struct powgov_runtime *runtime, struct phase_profile *updated_phase)
{
	// we have seen this before but we don't know if control flow will be the same
	// create a preliminary for this
	struct l3_graph_node *prelim;
	char rejected = 0;
	if (runtime->sampler->l3->prelim != NULL) // there is an existing prelim
	{
		// if exists itr->phase == prelim->phase && itr->next->phase == updated phase
		// 		we have seen this before, update pointers
		struct l3_graph_node *match = find_graph_node(runtime, updated_phase);
		if (match != NULL)
		{
			printf("reattach node for phase %d\n", updated_phase - runtime->classifier->phases);
			// reject the current preliminary since we have seen this before
			reject_preliminary(runtime, match);
			// current_node was pointing to the preliminary, but is now invalid
			rejected = 1;
			// then create a new preliminary for the updated phase
		}
		// else
		// 		we have not seen this before, commit prelim and create a new one
			//  we don't need to update pointers since the preliminary was not rejected
			//  fall through and update prelim->next with current_node pointer
	}
	// else
		// there is no prelim so a new node must have just been added
		// so we can just create a new prelim 
	// for all of the above, just create a new prelim
	prelim = get_new_graph_node(runtime);
	prelim->phase = updated_phase;
	prelim->next = NULL;
	runtime->sampler->l3->prelim = prelim;

	// advance current_node pointer
	if (rejected == 0)
	{
		runtime->sampler->l3->current_node->next = prelim;
	}
	// else
		// we reattached, so this edge already exists
	runtime->sampler->l3->current_node = prelim;

	dump_graph(runtime);
}
