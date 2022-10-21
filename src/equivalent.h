#ifndef AMULET2_SRC_EQUIVALENT_H_
#define AMULET2_SRC_EQUIVALENT_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"

#include <gmp.h>
#include <map>
#include <cmath>
#include <fstream>
#include <set>

typedef std::pair<unsigned int, int> pair_equiv;
static std::map<std::string, std::vector<int>> gateSim;
static std::map<std::string, pair_equiv> equivClass;
static std::set<std::string> node_name;
static int max_depth;
extern std::map<std::string, pair_equiv> get_equivClass;

void simulate();
void build_simulate_vector();
int check_equiv_candidate(std::string, std::string);

// build equivalence class
void SBIF(int windowDepth);
int valid_equiv(Gate *g, Gate *b, int windowDepth, int flag);
void valid_equiv_boolector(Gate *g, Gate *b, int windowDepth, int flag);

void gen_node_declare(Gate *g, int windowDepth);
void gen_node_formula(Gate *g, int windowDepth);

std::vector<Polynomial *> gen_equiv_poly(int window);

void print_dependency_all();
void print_dependency_gate(Gate *g, int depth);

#endif // AMULET2_SRC_EQUIVALENT_H_