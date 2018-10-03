#pragma once
#include "powgov.h"

#define MSR_TURBO_RATIO_LIMIT 0x1AD 

void dump_rapl(FILE *out);
void enable_turbo(struct powgov_runtime *runtime);
void disable_turbo(struct powgov_runtime *runtime);
void set_turbo_limit(unsigned int limit);
void dump_rapl_info(double power_unit);
void activate_performance_counters(struct powgov_runtime * runtime);
void set_perf(struct powgov_runtime *runtime, const unsigned freq);
void set_rapl(double sec, double watts, double pu, double su, unsigned affinity);
void set_rapl2(unsigned sec, double watts, double pu, double su, unsigned affinity);
void hwpstuff(FILE *out);
