

#include "SpliceMap_utils.h"

struct out_group_t {
	uint_fast32_t start;
	uint_fast32_t end;
	float val;
};


void write_coverage(ofstream &coverage, out_group_t &g, float curr_val,uint_fast32_t curr_loc, string chr);
void print_usage_and_exit();

