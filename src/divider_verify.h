#ifndef AMULET2_SRC_DIVIDER_VERIFY_H_
#define AMULET2_SRC_DIVIDER_VERIFY_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
#include "equivalent.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness_divider;
/*------------------------------------------------------------------------*/

static std::map<std::string, unsigned int> gateMap;

void dividerVerify(int split_num,
                   int window,
                   const char *inp_f = 0,
                   const char *out_f1 = 0,
                   const char *out_f2 = 0);

void split_solver(int split_num, int window, mpz_t low, mpz_t upper, Polynomial *spec, std::map<std::string, int> fix_q_value);

void dividerVerify_backup(int window,
                          const char *inp_f = 0,
                          const char *out_f1 = 0,
                          const char *out_f2 = 0);

void print_circuit_poly_divider(FILE *file);

// defines the divider specification
// \sum 2^(NN/2 - 1 - i) * bi =
//     (\sum 2^(NN/2 - 1 - i) * ai) * (\sum 2^(NN/2 - 1 - i) * si)
//     + \sum 2^(NN/2 - 1 - i) * s(i + NN/2)
Polynomial *print_spec_poly_divider(FILE *file);

// if lt(p) / lt(pi) return negative quotien
// else return 0
Polynomial *divide_by_lt(const Polynomial *p, const Monomial *tm);

Polynomial *divide_by_lm(const Polynomial *p, const Monomial *tm);

const Polynomial *reduce_divider(FILE *file);

// (divider_order == 1)
// boolector order:
// input: in_op2[NN/2 - 1], in_op2[NN/2 - 2], ..., in_op2[1], in_op1[NN/2 - 1],
//     in_op2[0], in_op1[NN/2 -2], ..., in_op1[1], in_op1[0]
// output: Q[NN/2 - 1], ..., Q[0], R[NN / 2 - 1], ..., R[0]
void fix_order();

const std::vector<Polynomial *> gen_substitute_poly(int split_num, int window, mpz_t l, mpz_t m, std::map<std::string, int> fix_q_value);
Polynomial *reduce_base(Polynomial *p, std::vector<Polynomial *> poly_set);
Polynomial *reduce_base_back(Polynomial *p, std::vector<Polynomial *> poly_set);

Polynomial *reduce_with_valid(Polynomial *p, std::vector<Polynomial *> poly_set, mpz_t l, mpz_t u);
const std::vector<Polynomial *> gen_substitute_poly_candidate(int split_num, int window, mpz_t l, mpz_t m, std::map<std::string, int> fix_q_value);

bool validate_rem(Polynomial *rem);

// Term *remainder_term(Term *t, Term *tv);
void remove_internal_xor_gates_divider(FILE *file);

bool same_poly(Polynomial *a, Polynomial *b);

void construct_gate_map();
const std::vector<Polynomial *> gen_original_poly();

#endif // AMULET2_SRC_DIVIDER_VERIFY_H_