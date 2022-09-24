#ifndef AMULET2_SRC_DIVIDER_VERIFY_H_
#define AMULET2_SRC_DIVIDER_VERIFY_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness_divider;
/*------------------------------------------------------------------------*/

void dividerVerify(const char *inp_f = 0,
                   const char *out_f1 = 0,
                   const char *out_f2 = 0);

void print_circuit_poly(FILE *file);

// defines the divider specification
// \sum 2^i * ai = (\sum 2^i * bi) * (\sum 2^i * si) + \sum 2^i * si+NN/2
Polynomial *print_spec_poly_divider(FILE *file);

const Polynomial *reduce_divider(FILE *file);

#endif // AMULET2_SRC_DIVIDER_VERIFY_H_