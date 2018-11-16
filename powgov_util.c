#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "msr_core.h"
#include "cpuid.h"
#include "msr_counters.h"
#include "master.h"
#include "powgov_util.h"

void dump_rapl(FILE *out)
{
	uint64_t rapl;
	read_msr_by_coord(0, 0, 0, MSR_PKG_POWER_LIMIT, &rapl);
	fprintf(out, "\trapl is %lx\n", (unsigned long) rapl);
}

void enable_turbo(struct powgov_runtime *runtime)
{
	uint64_t perf_ctl;
	read_msr_by_coord(0, 0, 0, IA32_PERF_CTL, &perf_ctl);
	perf_ctl &= 0xFFFFFFFEFFFFFFFFul;
	int i;
	for (i = 0; i < runtime->sys->num_cpu; i++)
	{
		write_msr_by_coord(0, i, 0, IA32_PERF_CTL, perf_ctl);
		/* write_msr_by_coord(1, i, 0, IA32_PERF_CTL, perf_ctl); */
	}
}

void disable_turbo(struct powgov_runtime *runtime)
{
	uint64_t perf_ctl;
	read_msr_by_coord(0, 0, 0, IA32_PERF_CTL, &perf_ctl);
	perf_ctl |= 0x0000000100000000ul;
	int i;
	for (i = 0; i < runtime->sys->num_cpu; i++)
	{
		write_msr_by_coord(0, i, 0, IA32_PERF_CTL, perf_ctl);
		/* write_msr_by_coord(1, i, 0, IA32_PERF_CTL, perf_ctl); */
	}
}

void set_turbo_limit(unsigned int limit)
{
	uint64_t turbo_limit;
	limit &= 0xFF;
	turbo_limit = 0x0 | (limit) | (limit << 8) | (limit << 16) | (limit << 24);
	//printf("set turbo limit %lx\n", (unsigned long) turbo_limit);
	write_msr_by_coord(0, 0, 0, MSR_TURBO_RATIO_LIMIT, turbo_limit);
}

void dump_rapl_info(double power_unit)
{
	uint64_t rapl_info;
	read_msr_by_coord(0, 0, 0, MSR_PKG_POWER_INFO, &rapl_info);
	fprintf(stderr, "\tRAPL INFO:\n\tTDP %lf (raw %lx)\n",
		(rapl_info & 0xEF) * power_unit, rapl_info & 0xEF);
}

void activate_performance_counters(struct powgov_runtime * runtime)
{
	set_all_pmc_ctrl(0x0, 0x43, 0x41, 0x2E, 1); // LLC miss
	set_all_pmc_ctrl(0x0, 0x43, 0x01, 0xA2, 2); // resource stalls
	set_all_pmc_ctrl(0x0, 0x43, 0x04, 0xA3, 3); // execution stalls
	set_all_pmc_ctrl(0x0, 0x43, 0x00, 0xC4, 4); // branch instructions retired
	//set_all_pmc_ctrl(0x0, 0x43, 0x02, 0xC7, 3); // SSE/AVX single precision retired
	//set_all_pmc_ctrl(0x0, 0x43, 0x01, 0xC7, 4); // SSE/AVX double precision retired
	//set_all_pmc_ctrl(0x0, 0x43, 0x04, 0xC5, 3); // branch misses retired
	enable_pmc();
	// enable fixed counters
	uint64_t ovf_ctrl;
	int ctr;
	for (ctr = 0; ctr < runtime->sys->num_cpu; ctr++)
	{
		uint64_t perf_global_ctrl;
		uint64_t fixed_ctr_ctrl;
		read_msr_by_coord(0, ctr, 0, IA32_PERF_GLOBAL_CTRL, &perf_global_ctrl);
		read_msr_by_coord(0, ctr, 0, IA32_FIXED_CTR_CTRL, &fixed_ctr_ctrl);
		write_msr_by_coord(0, ctr, 0, IA32_PERF_GLOBAL_CTRL, perf_global_ctrl | (0x7ul << 32) | 0x3);
		write_msr_by_coord(0, ctr, 0, IA32_FIXED_CTR_CTRL, fixed_ctr_ctrl | (0x3));
		write_msr_by_coord(0, ctr, 0, IA32_FIXED_CTR0, 0x0ul);
		read_msr_by_coord(0, ctr, 0, IA32_PERF_GLOBAL_OVF_CTRL, &ovf_ctrl);
		write_msr_by_coord(0, ctr, 0, IA32_PERF_GLOBAL_OVF_CTRL, ovf_ctrl & 0xFFFFFFFFFFFFFFFE);
		//write_msr_by_coord(0, ctr, 0, MSR_PKG_PERF_STATUS, 0);
	}
}

// this sets the processor p-state (frequency)
void set_perf(struct powgov_runtime *runtime, const unsigned freq)
{
	/* printf("received freq %u to set\n", freq); */
	uint64_t perf_ctl = 0x0ul;
	uint64_t freq_mask = freq;
	read_msr_by_coord(0, 0, 0, IA32_PERF_CTL, &perf_ctl);
	perf_ctl &= 0xFFFFFFFFFFFF0000ul;
	freq_mask <<= 8;
	perf_ctl |= freq_mask;
	/* printf("setting %lx\n", perf_ctl); */
	//write_msr_by_coord(0, tid, 0, IA32_PERF_CTL, perf_ctl);
	//write_msr_by_coord(0, tid, 1, IA32_PERF_CTL, perf_ctl);
	int i;
	for (i = 0; i < runtime->sys->num_cpu; i++)
	{
		write_msr_by_coord(0, i, 0, IA32_PERF_CTL, perf_ctl);
		//write_msr_by_coord(0, i, 1, IA32_PERF_CTL, perf_ctl);
		// TODO: should do this if second socket exists
		//write_msr_by_coord(1, i, 0, IA32_PERF_CTL, perf_ctl);
		//write_msr_by_coord(1, i, 1, IA32_PERF_CTL, perf_ctl);
	}
}

void set_rapl(double sec, double watts, double pu, double su, unsigned affinity)
{
	uint64_t power = (unsigned long) (watts / pu);
	uint64_t seconds;
	uint64_t timeval_y = 0, timeval_x = 0;
	double logremainder = 0;

	printf("setting RAPL to %lfwatts, %lf seconds\n", watts, sec);

	timeval_y = (uint64_t) log2(sec / su);
	// store the mantissa of the log2
	logremainder = (double) log2(sec / su) - (double) timeval_y;
	timeval_x = 0;
	// based on the mantissa, we can choose the appropriate multiplier
	if (logremainder > 0.15 && logremainder <= 0.45)
	{
		timeval_x = 1;
	}
	else if (logremainder > 0.45 && logremainder <= 0.7)
	{
		timeval_x = 2;
	}
	else if (logremainder > 0.7)
	{
		timeval_x = 3;
	}
	// store the bits in the Intel specified format
	seconds = (uint64_t) (timeval_y | (timeval_x << 5));
	uint64_t rapl = 0x0 | power | (seconds << 17);
	uint64_t oldrapl = 0x0;
	read_msr_by_coord(0, 0, 0, MSR_PKG_POWER_LIMIT, &oldrapl);

	rapl |= (1LL << 15) | (1LL << 16) | (oldrapl & 0xFFFFFFFF00000000);
	write_msr_by_coord(0, 0, 0, MSR_PKG_POWER_LIMIT, rapl);
}

void set_rapl2(unsigned sec, double watts, double pu, double su, unsigned affinity)
{
	uint64_t power = (unsigned long) (watts / pu);
	uint64_t seconds;

	seconds = sec;
	uint64_t rapl = 0x0 | power | (seconds << 17);
	rapl |= (1LL << 15) | (1LL << 16);

	uint64_t oldrapl;
	read_msr_by_coord(0, 0, 0, MSR_PKG_POWER_LIMIT, &oldrapl);

	rapl = (rapl << 32) | (oldrapl & 0x00000000ffffffff);
	write_msr_by_coord(0, 0, 0, MSR_PKG_POWER_LIMIT, rapl);
}

// this will check if HWP is enabled, power-governor is not compatible with HWP
void hwpstuff(FILE *out)
{
	uint64_t rax, rbx, rcx, rdx;
	rax = rbx = rcx = rdx = 0;
	cpuid(0x6, &rax, &rbx, &rcx, &rdx);
	if ((rax & (0x1ul << 7)) == 0)
	{
		fprintf(out, "\thwp is not supported or disabled\n");
		return;
	}

	uint64_t hwp = 0x0;
	read_msr_by_coord(0, 0, 0, 0x770, &hwp);
	fprintf(out, "\thwp %lx %s\n", hwp, (hwp == 0 ? "disabled" : "enabled"));
	if (hwp != 0)
	{
		uint64_t hwp_cap = 0x0;
		read_msr_by_coord(0, 0, 0, 0x771, &hwp_cap);
		fprintf( out, "\thwp cap %lx\n", hwp_cap);
		uint64_t hwp_req = 0x0;
		read_msr_by_coord(0, 0, 0, 0x772, &hwp_req);
		fprintf( out, "\thwp req %lx\n", hwp_req);
		uint64_t hwp_int = 0x0;
		read_msr_by_coord(0, 0, 0, 0x773, &hwp_int);
		fprintf( out, "\thwp int %lx\n", hwp_int);
		uint64_t hwp_log_req = 0x0;
		read_msr_by_coord(0, 0, 0, 0x774, &hwp_log_req);
		fprintf( out, "\thwp log req %lx\n", hwp_log_req);
		uint64_t hwp_stat = 0x0;
		read_msr_by_coord(0, 0, 0, 0x777, &hwp_stat);
		fprintf( out, "\thwp stat %lx\n", hwp_stat);

		hwp_req = 0x14150Alu;
		hwp_log_req = 0x14150Alu;
		//write_msr_by_coord(0, 0, 0, 0x772, hwp_req);
		//write_msr_by_coord(0, 0, 0, 0x772, hwp_req);
		read_msr_by_coord(0, 0, 0, 0x772, &hwp_req);
		fprintf( out, "\thwp req %lx\n", hwp_req);
		read_msr_by_coord(0, 0, 0, 0x774, &hwp_log_req);
		fprintf( out, "\thwp log req %lx\n", hwp_log_req);
	}
}
