#include "powgov_l2.h"
#include "powgov_l1.h"

void l2_analysis(struct powgov_runtime *runtime)
{
	runtime->sampler->l2->prev_sample = runtime->sampler->l2->new_sample;
	runtime->sampler->l2->new_sample = runtime->sampler->l1->new_sample;
	// TODO: software turbo
	// TODO: overpower
	
	double window_time = (runtime->sampler->l2->new_sample.tsc_data -
			runtime->sampler->l2->prev_sample.tsc_data) /
			((((runtime->sampler->l2->new_sample.frq_data & 0xFFFFul) >> 8) / 10.0) *
			 1000000000.0);
	double instant_time = (runtime->sampler->l1->new_sample.tsc_data
			- runtime->sampler->l1->prev_sample.tsc_data) /
			((((runtime->sampler->l1->new_sample.frq_data & 0xFFFFul) >> 8) / 10.0) *
			 1000000000.0);
	unsigned long window_diff = (runtime->sampler->l2->new_sample.energy_data -
		runtime->sampler->l2->prev_sample.energy_data);
	unsigned long instant_diff = (runtime->sampler->l1->new_sample.energy_data -
		runtime->sampler->l1->prev_sample.energy_data);
	double window_pow = window_diff * runtime->sys->rapl_energy_unit / window_time;
	double instant_pow = instant_diff * runtime->sys->rapl_energy_unit / instant_time;

	if (window_pow >= runtime->power->rapl1 || window_pow >= runtime->power->rapl2)
	{
		// we are up against the power limit, or over it. We should throttle more.
		if (instant_pow >= runtime->power->rapl1 || instant_pow >= runtime->power->rapl2)
		{
			printf("power excursion\n");
			runtime->power->excursion = 1;
		}
	}
}
