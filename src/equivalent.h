#ifndef AMULET2_SRC_EQUIVALENT_H_
#define AMULET2_SRC_EQUIVALENT_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
#include "parser.h"

#include <gmp.h>
#include <map>
#include <cmath>
#include <fstream>
#include <set>

typedef std::pair<unsigned int, int> pair_equiv;
static std::map<std::string, std::vector<int>> gateSim;
static std::map<std::string, int> tmpSim;
static std::map<std::string, pair_equiv> equivClass;
static std::set<std::string> node_name;
static int max_depth;
extern std::map<std::string, pair_equiv> get_equivClass;

void simulate_input();

// R^(0) < D * 2^(n-1)
bool constraint_input(mpz_t l, mpz_t u);
std::string concat_str_divident();
std::string concat_str_divider();
void simulate(mpz_t l, mpz_t u);
void build_simulate_vector(mpz_t l, mpz_t u);
int check_equiv_candidate(std::string, std::string);

// build equivalence class
void SBIF(int windowDepth, mpz_t l, mpz_t u);
int valid_equiv(Gate *g, Gate *b, int windowDepth, int flag, mpz_t l, mpz_t u);
void valid_equiv_boolector(Gate *g, Gate *b, int windowDepth, int flag);

void gen_node_declare(Gate *g, int windowDepth);
void gen_node_formula(Gate *g, int windowDepth);

std::vector<Polynomial *> gen_equiv_poly(int window, mpz_t l, mpz_t m);

void print_dependency_all();
void print_dependency_gate(Gate *g, int depth);

#endif // AMULET2_SRC_EQUIVALENT_H_