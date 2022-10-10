#include "equivalent.h"

void simulate()
{
    // std::map<std::string, int> gate_value;
    for (unsigned i = 0; i < NN; i++)
    {
        const Var *x = gates[i]->get_var();
        int assign = rand() % 2;
        gateSim[x->get_name()].emplace_back(assign);
    }

    for (unsigned i = NN; i < num_gates; i++)
    {
        Gate *g = gates[i];
        std::string g_name = g->get_var_name();
        // std::cout << g->get_var_name() << " ";
        if (i < M - 1)
        {
            aiger_and *and1 = is_model_and(g->get_var_num());
            unsigned l = and1->rhs0, r = and1->rhs1;
            Gate *l_gate = gate(l), *r_gate = gate(r);
            std::string l_name = l_gate->get_var_name();
            std::string r_name = r_gate->get_var_name();
            int l_value, r_value;
            // std::cout << aiger_sign(l) << " " << l_gate->get_var_name()
            //           << " " << aiger_sign(r) << " " << r_gate->get_var_name() << '\n';

            l_value = (aiger_sign(l) == 1) ? 1 - gateSim[l_name].back() : gateSim[l_name].back();
            r_value = (aiger_sign(r) == 1) ? 1 - gateSim[r_name].back() : gateSim[r_name].back();
            gateSim[g_name].emplace_back(l_value * r_value);
        }
        else
        {
            assert(g->get_output());

            Gate *output_child_gate = g->children_front();
            std::string output_child = output_child_gate->get_var_name();
            int output_child_value = (aiger_sign(slit(i - M + 1)) == 1) ? 1 - gateSim[output_child].back() : gateSim[output_child].back();
            gateSim[g_name].emplace_back(output_child_value);
            // std::cout << aiger_sign(slit(i - M + 1)) << " " << model_output_gate->get_var_name() << "\n";
        }
    }
}

void build_simulate_vector()
{
    unsigned int simulate_num = pow(2, NN) / 2;
    srand(time(NULL));
    for (unsigned i = 0; i < simulate_num; i++)
    {
        simulate();
    }
}

int check_equiv_candidate(std::string gate1_name, std::string gate2_name)
{
    int equiv = 1;
    int antiEquiv = 1;
    std::vector<int> gate_vector1 = gateSim[gate1_name];
    std::vector<int> gate_vector2 = gateSim[gate2_name];
    int vecSize = gate_vector1.size();
    // for (int i = 0; i < vecSize; i++)
    // {
    //     std::cout << gate_vector1[i] << " ";
    // }
    // std::cout << "\n";
    // vecSize = gate_vector2.size();
    // for (int i = 0; i < vecSize; i++)
    // {
    //     std::cout << gate_vector2[i] << " ";
    // }
    // std::cout << "\n";
    for (int i = 0; i < vecSize; i++)
    {
        if (gate_vector1[i] == gate_vector2[i])
            antiEquiv = 0;
        if (gate_vector1[i] != gate_vector2[i])
            equiv = 0;
    }

    if (equiv == 1)
        return 1;
    if (antiEquiv == 1)
        return -1;
    else
        return 0;
}

void SBIF(int windowDepth)
{
    // initial
    build_simulate_vector();

    for (unsigned i = NN; i < num_gates; i++)
    {
        std::string gate_name = gates[i]->get_var_name();
        equivClass[gate_name] = pair_equiv(i, true);
    }

    // update equivalence class
    for (unsigned i = NN; i < M - 1; i++)
    {
        std::string a_name = gates[i]->get_var_name();

        for (unsigned j = NN; j < i; j++)
        {
            std::string b_name = gates[j]->get_var_name();
            int equiv = check_equiv_candidate(a_name, b_name);
            // std::cout << a_name << " " << b_name << " equiv: " << equiv << "\n";
            if (equiv == 0)
                continue;
            assert(equiv == 1 || equiv == -1);

            // std::cout << "equivClass[" << a_name << "]: " << equivClass[a_name] << " equivClass[" << b_name << "]: "
            //           << equivClass[b_name] << "\n";
            // if a \notin [b]
            if (equivClass[a_name].first != equivClass[b_name].first)
            {
                valid_equiv(gates[i], gates[j], windowDepth, equiv);
                // std::cout << a_name << " " << b_name << " : " << equiv_result << "\n";
            }
        }
    }
    // for (unsigned i = NN; i < M - 1; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     std::cout << equivClass[a_name].first << "\n";
    // }
}

void valid_equiv(Gate *a, Gate *b, int windowDepth, int flag)
{
    std::ofstream fp;
    fp.open("valid.smt2");
    fp << "(set-logic ALL)\n";
    node_name.clear();

    gen_node_declare(a, windowDepth);
    gen_node_declare(b, windowDepth);

    std::string a_name = a->get_var_name();
    std::string b_name = b->get_var_name();
    node_name.insert(a_name);
    node_name.insert(b_name);
    for (auto n : node_name)
    {
        fp << "(declare-fun " << n << "() Bool)\n";
    }
    fp.close();
    gen_node_formula(a, windowDepth);
    gen_node_formula(b, windowDepth);

    fp.open("valid.smt2", std::ios::app);

    if (flag == 1)
    {
        // a xor b = ((not a) and (b)) or ((a) and (not b))
        // if sat, there exists an assignment that prove a and b not equiv/anti-equiv
        fp << "(assert (or (and (not " << a_name << ") " << b_name << ") (and " << a_name << " (not " << b_name << "))))\n";
        fp << "(check-sat)\n";
    }

    else
    {
        assert(flag == -1);
        fp << "(assert (or (and (not " << a_name << ") "
           << "(not " << b_name << ")) (and " << a_name << " " << b_name << ")))\n";
        fp << "(check-sat)\n";
    }

    fp.close();

    std::string command = "boolector valid.smt2 > valid_output";
    system(command.c_str());

    std::ifstream FILE("valid_output");
    std::string line;
    getline(FILE, line);
    if (line == "unsat")
    {
        std::cout << a_name << " and " << b_name << " are (anti-)equivlent\n";
        if (equivClass[a_name].first < equivClass[b_name].first)
        {
            if (flag == 1)
            {
                equivClass[b_name] = pair_equiv(equivClass[a_name].first, true);
            }
            else if (flag == -1)
            {
                equivClass[b_name] = pair_equiv(equivClass[a_name].first, false);
            }
        }
        else
        {
            if (flag == 1)
            {
                equivClass[a_name] = pair_equiv(equivClass[b_name].first, true);
            }
            else if (flag == -1)
            {
                equivClass[a_name] = pair_equiv(equivClass[b_name].first, false);
            }
        }
        // return 1;
    }
    else if (line != "sat")
    {
        msg("Unexpect SMT solver output!");
    }

    // return 0;
}

void gen_node_declare(Gate *g, int windowDepth)
{
    if (windowDepth == 0)
        return;
    if (g->get_input())
        return;
    else
    {
        aiger_and *and1 = is_model_and(g->get_var_num());
        unsigned l = and1->rhs0, r = and1->rhs1;
        Gate *l_gate = gate(l), *r_gate = gate(r);
        // std::cout << "l: " << l_gate->get_var_name() << "\n";
        // std::cout << "r: " << r_gate->get_var_name() << "\n";
        node_name.insert(l_gate->get_var_name());
        node_name.insert(r_gate->get_var_name());
        gen_node_declare(l_gate, windowDepth - 1);
        gen_node_declare(r_gate, windowDepth - 1);
    }
}

void gen_node_formula(Gate *g, int windowDepth)
{
    if (windowDepth == 0)
        return;
    if (g->get_input())
        return;
    else
    {
        std::string g_name = g->get_var_name();
        aiger_and *and1 = is_model_and(g->get_var_num());
        unsigned l = and1->rhs0, r = and1->rhs1;
        Gate *l_gate = gate(l), *r_gate = gate(r);
        std::string l_name = l_gate->get_var_name();
        std::string r_name = r_gate->get_var_name();
        std::string l_literal;
        std::string r_literal;

        /*
         * ; Variable declarations
         *(declare-fun a () Bool)
         *(declare-fun b () Bool)
         *
         *; Constraints
         *(assert (and (not a) (not b)))
         *(assert (and a (not b)))
         *
         *; Solve
         *(check-sat)
         *(get-model)
         */

        if (aiger_sign(l) == 0)
            l_literal = l_name;
        else
            l_literal = "(not " + l_name + ")";
        if (aiger_sign(r) == 0)
            r_literal = r_name;
        else
            r_literal = "(not " + r_name + ")";
        // fprintf(f, "(assert (and %s %s))\n", l_literal.c_str(), r_literal.c_str());

        std::ofstream fp;
        fp.open("valid.smt2", std::ios::app);
        fp << "(assert (= " << g_name << " (and " << l_literal << " " << r_literal << ")))\n";
        fp.close();

        gen_node_formula(l_gate, windowDepth - 1);
        gen_node_formula(r_gate, windowDepth - 1);
    }
}

std::vector<Polynomial *> gen_equiv_poly(int window)
{
    SBIF(window);

    std::vector<Polynomial *> poly_set;

    for (unsigned i = NN; i < num_gates; i++)
    {
        Gate *g = gates[i];
        pair_equiv equiv_tag = equivClass[g->get_var_name()];
        // std::cout << equiv_tag.first << "\n";
        // std::cout << equiv_tag.second << "\n";
        if (equiv_tag.first != i)
        {
            Gate *gs = gates[equiv_tag.first];
            // std::cout << "Should substitute " << g->get_var_name() << " to " << gs->get_var_name() << "\n";
            if (equiv_tag.second == true)
            {
                const Var *v1 = g->get_var();
                Term *t1 = new_term(v1, 0);
                Monomial *m1 = new Monomial(minus_one, t1);
                push_mstack_end(m1);
                const Var *v2 = gs->get_var();
                Term *t2 = new_term(v2, 0);
                Monomial *m2 = new Monomial(one, t2);
                push_mstack_end(m2);
            }
            else
            {
                assert(equiv_tag.second == false);
                const Var *v1 = g->get_var();
                Term *t1 = new_term(v1, 0);
                Monomial *m1 = new Monomial(minus_one, t1);
                push_mstack_end(m1);
                const Var *v2 = gs->get_var();
                Term *t2 = new_term(v2, 0);
                Monomial *m2 = new Monomial(minus_one, t2);
                push_mstack_end(m2);
                Monomial *constant_one = new Monomial(one, 0);
                push_mstack_end(constant_one);
            }

            Polynomial *p_equiv = build_poly();
            poly_set.emplace_back(p_equiv);
        }
    }

    return poly_set;
}