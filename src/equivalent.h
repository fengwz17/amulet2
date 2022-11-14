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
#include <z3++.h>

typedef std::pair<unsigned int, int> pair_equiv;
typedef std::pair<std::string, z3::expr> pair_expr;
static std::map<std::string, std::vector<int>> gateSim;
static std::map<std::string, int> tmpSim;
static std::map<std::string, pair_equiv> equivClass;
static std::set<std::string> node_name;
static int max_depth;
extern std::map<std::string, pair_equiv> get_equivClass;
static int same;
// std::vector<std::vector<int>> input_vec;
static std::map<std::string, z3::expr> name_expr;

void simulate_input();
void simulate_input_vector(int split_num, mpz_t l, mpz_t u);

// R^(0) < D * 2^(n-1)
bool constraint_input(mpz_t l, mpz_t u);
std::string concat_str_divident();
std::string concat_str_divider();
void simulate(int split_num, mpz_t l, mpz_t u);
void build_simulate_vector(int split_num, mpz_t l, mpz_t u);
int check_equiv_candidate(std::string, std::string);

// build equivalence class
void SBIF(int split_num, int windowDepth, mpz_t l, mpz_t u);
void SBIF_candidate(int split_num, int windowDepth, mpz_t l, mpz_t u);

int valid_equiv(Gate *g, Gate *b, int windowDepth, int flag, mpz_t l, mpz_t u);
// void valid_equiv_boolector(Gate *g, Gate *b, int windowDepth, int flag);
int valid_equiv_z3(Gate *a, Gate *b, int windowDepth, int flag, mpz_t l, mpz_t u, z3::context *ctx);

void gen_node_declare(Gate *g, int windowDepth);
void gen_node_formula(Gate *g, int windowDepth);
void gen_node_formula_z3(Gate *g, int windowDepth, z3::solver *solver, z3::context *ctx);

std::vector<Polynomial *> gen_equiv_poly(int split_num, int window, mpz_t l, mpz_t m);
std::vector<Polynomial *> gen_equiv_poly_candidate(int split_num, int window, mpz_t l, mpz_t m);

void print_dependency_all();
void print_dependency_gate(Gate *g, int depth);

void set_equiv(std::string aName, std::string bName, int equiv_flag);

#endif // AMULET2_SRC_EQUIVALENT_H_