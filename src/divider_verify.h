#ifndef AMULET2_SRC_DIVIDER_VERIFY_H_
#define AMULET2_SRC_DIVIDER_VERIFY_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
#include "parser.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness_divider;
/*------------------------------------------------------------------------*/

void dividerVerify(const char *inp_f = 0,
                   const char *out_f1 = 0,
                   const char *out_f2 = 0);

void print_circuit_poly_divider(FILE *file);

// defines the divider specification
// \sum 2^i * ai = (\sum 2^i * bi) * (\sum 2^i * si) + \sum 2^i * si+NN/2
Polynomial *print_spec_poly_divider(FILE *file);

// if lt(p) / lt(pi) return negative quotien
// else return 0
Polynomial *divide_by_lt(const Polynomial *p, const Term *t);

const Polynomial *reduce_divider(FILE *file);

void fix_order();

#endif // AMULET2_SRC_DIVIDER_VERIFY_H_