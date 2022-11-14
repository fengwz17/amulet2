#include "equivalent.h"

void simulate(int split_num, mpz_t l, mpz_t u)
{

    // simulate_input();
    // if (divider_order == 2)
    // {
    //     while (!constraint_input(l, u))
    //     {
    //         simulate_input();
    //     }
    // }

    // sample an input that satisfies the constrint
    simulate_input_vector(split_num, l, u);

    std::vector<int> tmp_input;
    tmp_input.clear();
    for (unsigned i = 0; i < NN; i++)
    {
        const Var *x = gates[i]->get_var();
        gateSim[x->get_name()].emplace_back(tmpSim[x->get_name()]);

        tmp_input.emplace_back(tmpSim[x->get_name()]);
    }
    // input_vec.emplace_back(tmp_input);

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

void simulate_input()
{
    for (unsigned i = 0; i < NN; i++)
    {
        const Var *x = gates[i]->get_var();
        int assign = rand() % 2;
        tmpSim[x->get_name()] = assign;
    }

    // for (unsigned i = 0; i < NN; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     int av = tmpSim[a_name];

    //     std::cout << av << " ";
    // }
    // std::cout << "\n";
}

void simulate_input_vector(int split_num, mpz_t l, mpz_t u)
{
    clock_t seed = clock();
    // unsigned int seed = time(NULL);
    mpz_t D;
    mpz_t R0;
    mpz_init(D);
    mpz_init(R0);

    // D size is width - 1
    // sample D in range (0, 2^{width - 1})
    int width = NN / 3 + 1;
    char strD[width];
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, seed);

    // D != 0
    while (mpz_cmp_ui(D, 0) == 0)
        mpz_urandomb(D, rstate, width - 1);

    // gmp_printf("generate D: %Zd\n", D);

    mpz_get_str(strD, 2, D);
    int Dlen = strlen(strD);
    for (int i = width - 2; i >= 0; i--)
    {
        if (i <= Dlen - 1)
            tmpSim[gates[i]->get_var_name()] = strD[Dlen - 1 - i] - '0';
        else
            tmpSim[gates[i]->get_var_name()] = 0;
    }
    // for (int i = 0; i <= width - 2; i++)
    //     std::cout << tmpSim[gates[i]->get_var_name()] << " ";
    // std::cout << "\n";

    // l * D <= R0 < m * D
    // m - l = step = pow(2, width - split_num),
    // e.g., width = 4, split_num = 2, then range is [0, 2^2) and [2^2, 2^3), step = 2^2
    // sample R0 in range [0, step * D) + l
    char strR0[width];
    mpz_t step;
    mpz_init(step);
    // R0 < step * D
    mpz_pow_ui(step, base, width - split_num);
    mpz_mul(step, step, D);
    mpz_urandomb(R0, rstate, width - split_num + Dlen);
    while (mpz_cmp(R0, step) != -1)
        mpz_urandomb(R0, rstate, width - split_num + Dlen);
    // gmp_printf("generate D: %Zd\n", D);
    // gmp_printf("generate step: %Zd\n", step);
    // gmp_printf("generate R0: %Zd\n", R0);

    mpz_add(R0, R0, l);
    mpz_get_str(strR0, 2, R0);
    int R0len = strlen(strR0);

    for (int i = 2 * width - 3; i >= 0; i--)
    {
        if (i <= R0len - 1)
            tmpSim[gates[i + width - 1]->get_var_name()] = strR0[R0len - 1 - i] - '0';
        else
            tmpSim[gates[i + width - 1]->get_var_name()] = 0;
    }
    // gmp_printf("generate R0: %Zd\n", R0);
    // for (int i = 0; i <= 2 * width - 3; i++)
    //     std::cout << tmpSim[gates[i + width - 1]->get_var_name()];
    // std::cout << "\n";

    gmp_randclear(rstate);
    mpz_clear(step);
    mpz_clear(D);
    mpz_clear(R0);
}

bool constraint_input(mpz_t l, mpz_t u)
{
    // for (unsigned i = 0; i < NN; i++)
    // {
    //     const Var *x = gates[i]->get_var();
    // }
    // int fm = 0;
    // int fl = 0;
    int width = NN / 3 + 1;

    mpz_t Dvalue;
    mpz_t R0value;
    mpz_t tmp_pow;
    mpz_init(Dvalue);
    mpz_init(R0value);
    mpz_init(tmp_pow);

    // D
    for (int i = width - 2; i >= 0; i--)
    {
        // std::cout << "tmpsim: " << tmpSim[gates[i]->get_var_name()] << "\n";
        if (tmpSim[gates[i]->get_var_name()] == 1)
        {
            mpz_pow_ui(tmp_pow, base, i);
            mpz_add(Dvalue, Dvalue, tmp_pow);
        }
        // gmp_printf("%Zd\n", Dvalue);
    }
    // gmp_printf("Dvalue: %Zd\n", Dvalue);

    // R0
    for (int i = 3 * width - 4; i >= width - 1; i--)
    {
        // std::cout << tmpSim[gates[i]->get_var_name()];
        if (tmpSim[gates[i]->get_var_name()] == 1)
        {
            mpz_pow_ui(tmp_pow, base, i - width + 1);
            mpz_add(R0value, R0value, tmp_pow);
        }
        // gmp_printf("%Zd\n", R0value);
    }
    // gmp_printf("R0value: %Zd\n", R0value);

    mpz_t lowerBound;
    mpz_t upperBound;
    mpz_init(lowerBound);
    mpz_init(upperBound);
    mpz_set(lowerBound, l);
    mpz_set(upperBound, u);

    // gmp_printf("l: %Zd\n", l);
    // gmp_printf("u: %Zd\n", u);
    // gmp_printf("lower: %Zd\n", lowerBound);
    // gmp_printf("upper: %Zd\n", upperBound);

    mpz_mul(lowerBound, Dvalue, lowerBound);
    mpz_mul(upperBound, Dvalue, upperBound);

    // gmp_printf("lower: %Zd\n", lowerBound);
    // gmp_printf("upper: %Zd\n", upperBound);

    int r = 0;

    // R0value < lowerBound return -1
    // = return 0
    // > return 1
    if (mpz_cmp(R0value, lowerBound) != -1 && mpz_cmp(R0value, upperBound) == -1)
        r = 1;

    mpz_clear(Dvalue);
    mpz_clear(R0value);
    mpz_clear(lowerBound);
    mpz_clear(upperBound);
    mpz_clear(tmp_pow);
    if (r == 1)
        return true;

    // std::cout << "false\n";
    return false;
}

void build_simulate_vector(int split_num, mpz_t l, mpz_t u)
{
    // srand(time(NULL));
    // unsigned int simulate_num = pow(2, NN) / 2;
    // int width = NN / 3 + 1;
    unsigned int simulate_num = pow(2, 10);
    // same = 0;
    // input_vec.clear();
    for (unsigned i = 0; i < simulate_num; i++)
    {
        simulate(split_num, l, u);
    }
    // for (int i = 0; i < simulate_num; i++)
    // {
    //     gmp_printf("upper: %Zd\n", D_vec[i]);
    //     gmp_printf("upper: %Zd\n", R_vec[i]);
    // }

    // for (int i = 0; i < input_vec.size(); i++)
    // {
    //     for (int j = 0; j < input_vec[i].size(); j++)
    //     {
    //         std::cout << input_vec[i][j];
    //     }
    //     std::cout << "\n";
    // }
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
        if (antiEquiv == 0 && equiv == 0)
            break;
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

void SBIF(int split_num, int windowDepth, mpz_t l, mpz_t u)
{
    double sim_time = process_time();

    // initial
    build_simulate_vector(split_num, l, u);

    // for (unsigned i = 0; i < NN; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     std::vector<int> avec = gateSim[a_name];
    //     for (auto i : avec)
    //     {
    //         std::cout << i << " ";
    //     }
    //     std::cout << "\n";
    // }

    double finish_sim_time = process_time();

    msg("used time for build_simulate_vec: %18.2f seconds",
        finish_sim_time - sim_time);

    for (unsigned i = 0; i < num_gates; i++)
    {
        std::string gate_name = gates[i]->get_var_name();
        equivClass[gate_name] = pair_equiv(i, 1);
    }

    double equi_time = process_time();

    // int count = 0;
    // update equivalence class
    int count_candidate = 0;
    int count_cannot_prove = 0;
    int each_candidate = 0;
    int count_signal = 0;
    int signal_flag = 0;
    int prove = 0;

    z3::context *ctx = new z3::context();
    for (unsigned i = NN; i < M - 1; i++)
    {
        std::string a_name = gates[i]->get_var_name();
        signal_flag = 0;
        each_candidate = 0;
        // for (unsigned j = NN; j < i; j++)
        for (unsigned j = i - 1; j >= NN; j--)
        {
            std::string b_name = gates[j]->get_var_name();
            int equiv = check_equiv_candidate(a_name, b_name);
            if (equiv == 0)
                continue;
            assert(equiv == 1 || equiv == -1);

            // std::cout << a_name << " " << b_name << " equiv: " << equiv << "\n";

            count_candidate++;
            each_candidate++;
            // std::cout << "equivClass[" << a_name << "]: " << equivClass[a_name].first << " " << equivClass[a_name].second
            //           << " equivClass[" << b_name << "]: "
            //           << equivClass[b_name].first << " " << equivClass[b_name].second << "\n";
            // // if a \notin [b]
            // std::cout << i << " " << j << " " << count_candidate << "\n";
            if (equivClass[a_name].first != equivClass[b_name].first)
            {
                // if (valid_equiv_z3(gates[i], gates[j], windowDepth, equiv, l, u, ctx) == 1)
                // if (valid_equiv(gates[i], gates[j], windowDepth, equiv, l, u) == 1)
                // {
                set_equiv(gates[i]->get_var_name(), gates[j]->get_var_name(), equiv);
                signal_flag = 1;
                prove++;
                break;
                // }
            }
            // // count++;
        }
        if (signal_flag == 0)
            count_signal++;
        if (signal_flag == 0)
            count_cannot_prove += each_candidate;
    }
    double end_equi_time = process_time();
    msg("used time for check_equiv: %18.2f seconds",
        end_equi_time - equi_time);
    delete (ctx);
    std::cout << count_candidate << "\n";
    std::cout << count_cannot_prove << "\n";
    // std::cout << count_signal << "\n";
    std::cout << prove << "\n";

    // std::cout << "find " << count << " candidates.\n";
    // for (unsigned i = NN; i < M - 1; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     std::cout << equivClass[a_name].first << "\n";
    // }
    // equivClassPoint = &equivClass;
    get_equivClass = equivClass;
}

void SBIF_candidate(int split_num, int windowDepth, mpz_t l, mpz_t u)
{
    double sim_time = process_time();

    // initial
    build_simulate_vector(split_num, l, u);

    // for (unsigned i = 0; i < NN; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     std::vector<int> avec = gateSim[a_name];
    //     for (auto i : avec)
    //     {
    //         std::cout << i << " ";
    //     }
    //     std::cout << "\n";
    // }

    double finish_sim_time = process_time();

    msg("used time for build_simulate_vec: %18.2f seconds",
        finish_sim_time - sim_time);

    for (unsigned i = 0; i < num_gates; i++)
    {
        std::string gate_name = gates[i]->get_var_name();
        equivClass[gate_name] = pair_equiv(i, 1);
    }

    double equi_time = process_time();

    // int count = 0;
    // update equivalence class
    int count_candidate = 0;
    int count_cannot_prove = 0;
    int each_candidate = 0;
    int count_signal = 0;
    int signal_flag = 0;
    int prove = 0;

    z3::context *ctx = new z3::context();
    for (unsigned i = NN; i < M - 1; i++)
    {
        std::string a_name = gates[i]->get_var_name();
        signal_flag = 0;
        each_candidate = 0;
        for (unsigned j = NN; j < i; j++)
        // for (unsigned j = i - 1; j >= NN; j--)
        {
            std::string b_name = gates[j]->get_var_name();
            int equiv = check_equiv_candidate(a_name, b_name);
            if (equiv == 0)
                continue;
            assert(equiv == 1 || equiv == -1);

            // std::cout << a_name << " " << b_name << " equiv: " << equiv << "\n";

            count_candidate++;
            each_candidate++;
            // std::cout << "equivClass[" << a_name << "]: " << equivClass[a_name].first << " " << equivClass[a_name].second
            //           << " equivClass[" << b_name << "]: "
            //           << equivClass[b_name].first << " " << equivClass[b_name].second << "\n";
            // // if a \notin [b]
            // std::cout << i << " " << j << " " << count_candidate << "\n";
            if (equivClass[a_name].first != equivClass[b_name].first)
            {
                // if (valid_equiv_z3(gates[i], gates[j], windowDepth, equiv, l, u, ctx) == 1)
                // if (valid_equiv(gates[i], gates[j], windowDepth, equiv, l, u) == 1)
                // {
                set_equiv(gates[i]->get_var_name(), gates[j]->get_var_name(), equiv);
                signal_flag = 1;
                // prove++;
                break;
                // }
            }
            // // count++;
        }
        if (signal_flag == 0)
            count_signal++;
        if (signal_flag == 0)
            count_cannot_prove += each_candidate;
    }
    double end_equi_time = process_time();
    msg("used time for check_equiv: %18.2f seconds",
        end_equi_time - equi_time);
    delete (ctx);
    std::cout << count_candidate << "\n";
    // std::cout << count_cannot_prove << "\n";
    // std::cout << count_signal << "\n";
    // std::cout << prove << "\n";

    // std::cout << "find " << count << " candidates.\n";
    // for (unsigned i = NN; i < M - 1; i++)
    // {
    //     std::string a_name = gates[i]->get_var_name();
    //     std::cout << equivClass[a_name].first << "\n";
    // }
    // equivClassPoint = &equivClass;
    get_equivClass = equivClass;
}

std::string concat_str_divident()
{
    int width = NN / 3 + 1;
    std::string concat_str = "(concat b1 b0)";
    // std::string concat_str = "(concat b" + std::to_string(width) + " b" + std::to_string(width - 1) + ")";

    for (int i = 2; i <= 2 * width - 3; i++)
    {
        concat_str = "(concat b" + std::to_string(i) + " " + concat_str + ")";
    }
    // for (int i = width + 1; i <= 2 * width - 3; i++)
    // {
    //     concat_str = "(concat b" + std::to_string(i) + " " + concat_str + ")";
    // }

    return concat_str;
}

std::string concat_str_divider()
{
    std::string concat_str = "(concat a1 a0)";
    int width = NN / 3 + 1;
    for (int i = 2; i <= width - 2; i++)
    {
        concat_str = "(concat a" + std::to_string(i) + " " + concat_str + ")";
    }
    std::string hi_bit = "#b";
    for (int i = 1; i <= width - 1; i++)
    {
        hi_bit = hi_bit + "0";
    }
    concat_str = "(concat " + hi_bit + " " + concat_str + ")";
    return concat_str;
}

int valid_equiv_z3(Gate *a, Gate *b, int windowDepth, int flag, mpz_t l, mpz_t u, z3::context *ctx)
{
    // std::ofstream fp;
    // fp.open("valid.smt2");
    // fp << "(set-logic QF_BV)\n";
    z3::solver *solver = new z3::solver(*ctx);
    // Z3_symbol logic_symbol = Z3_mk_string_symbol(*ctx, "QF_BV");
    // Z3_solver solver = Z3_mk_solver_for_logic(*ctx, logic_symbol);

    // z3::context ctx;
    // z3::solver solver(ctx);

    node_name.clear();
    int width = NN / 3 + 1;

    char *lc = mpz_get_str(NULL, 10, l);
    char *uc = mpz_get_str(NULL, 10, u);
    // char lc[2 * width - 2];
    // char uc[2 * width - 2];

    // for (int i = 0; i < 2 * width - 2; i++)
    // {
    //     if (i >= 2 * width - 2 - strlen(lc_tmp))
    //         lc[i] = lc_tmp[i - (2 * width - 2 - strlen(lc_tmp))];
    //     else
    //         lc[i] = '0';
    // }
    // for (int i = 0; i < 2 * width - 2; i++)
    // {
    //     if (i >= 2 * width - 2 - strlen(uc_tmp))
    //         uc[i] = uc_tmp[i - (2 * width - 2 - strlen(uc_tmp))];
    //     else
    //         uc[i] = '0';
    // }
    // for (int i = 0; i < 2 * width - 2; i++)
    // {
    //     std::cout << uc[i];
    // }
    // std::cout << "\n";

    if (divider_order == 2)
    {
        for (unsigned i = 0; i < NN; i++)
        {
            std::string input_name = gates[i]->get_var_name();
            node_name.insert(input_name);
        }
    }

    gen_node_declare(a, windowDepth);
    gen_node_declare(b, windowDepth);

    std::string a_name = a->get_var_name();
    std::string b_name = b->get_var_name();
    node_name.insert(a_name);
    node_name.insert(b_name);
    // for (auto n : node_name)
    // {
    //     fp << "(declare-fun " << n << "() (_ BitVec 1))\n";
    // }
    name_expr.clear();
    for (auto n : node_name)
    {
        z3::expr nameTemp = ctx->bv_const(n.c_str(), 1);
        // std::cout << nameTemp << " ";
        name_expr.insert(pair_expr(n, nameTemp));
    }

    /*
     * (set-logic ALL)
     * (declare-fun a0 () (_ BitVec 1))
     * (declare-fun a1 () (_ BitVec 1))
     * (declare-fun b0 () (_ BitVec 1))
     * (declare-fun b1 () (_ BitVec 1))
     * (declare-fun b2 () (_ BitVec 1))
     * (declare-fun Rn () (_ BitVec 3))
     * (declare-fun D () (_ BitVec 3))
     * (assert (= Rn (concat b2 (concat b1 b0))))
     * (assert (= D (concat (concat a1 a0) #b0)))
     * (assert (bvult Rn D))
     */
    if (divider_order == 2)
    {
        // std::string concat_str_R0;
        // concat_str_R0 = concat_str_divident();
        // std::string concat_str_D = concat_str_divider();
        // fp << "(declare-fun Rn () (_ BitVec " << 2 * width - 2 << "))\n";
        // fp << "(declare-fun Dl () (_ BitVec " << 2 * width - 2 << "))\n";
        // fp << "(declare-fun Du () (_ BitVec " << 2 * width - 2 << "))\n";
        // fp << "(assert (= Rn " << concat_str_R0 << "))\n";
        // fp << "(assert (= Dl "
        //    << "(bvmul " << concat_str_D << " (_ bv" << lc << " " << 2 * width - 2 << "))))\n";
        // fp << "(assert (= Du "
        //    << "(bvmul " << concat_str_D << " (_ bv" << uc << " " << 2 * width - 2 << "))))\n";
        // fp << "(assert (bvult Rn Du))\n";
        // fp << "(assert (bvuge Rn Dl))\n";

        z3::expr R0 = ctx->bv_const("R0", 2 * width - 2);
        z3::expr Dl = ctx->bv_const("Dl", 2 * width - 2);
        z3::expr Du = ctx->bv_const("Du", 2 * width - 2);
        z3::expr D = ctx->bv_const("D", 2 * width - 2);

        // extract(0,0) D = a0
        for (int i = 0; i < width - 1; i++)
        {
            z3::expr atmp = name_expr.at(gates[i]->get_var_name());
            solver->add(D.extract(i, i) == atmp);
        }

        // extract(0,0) R0 = b0
        for (int i = 0; i < 2 * width - 2; i++)
        {
            z3::expr btmp = name_expr.at(gates[i + width - 1]->get_var_name());
            solver->add(R0.extract(i, i) == btmp);
            // std::cout << name_expr[gates[i + width - 1]->get_var_name()];
        }

        z3::expr lexpr = ctx->bv_val(lc, 2 * width - 2);
        z3::expr uexpr = ctx->bv_val(uc, 2 * width - 2);
        // std::cout << "lexpr: " << lexpr << "\n";
        // std::cout << "uexpr: " << uexpr << "\n";
        solver->add(R0 < uexpr * D);
        solver->add(R0 >= lexpr * D);
    }

    // fp.close();

    // fp.open("valid.smt2", std::ios::app);

    z3::expr ac = name_expr.at(a_name);
    z3::expr bc = name_expr.at(b_name);

    if (flag == 1)
    {
        // a xor b = ((not a) and (b)) or ((a) and (not b))
        // if sat, there exists an assignment that prove a and b not equiv/anti-equiv
        // fp << "(assert (or (and (not " << a_name << ") " << b_name << ") (and " << a_name << " (not " << b_name << "))))\n";
        // fp << "(assert (= #b1 (bvxor " << a_name << " " << b_name << ")))\n";

        solver->add((ac ^ bc) == 1);
    }

    else
    {
        assert(flag == -1);
        // fp << "(assert (or (and (not " << a_name << ") "
        //    << "(not " << b_name << ")) (and " << a_name << " " << b_name << ")))\n";
        // fp << "(assert (= " << a_name << " " << b_name << "))\n";
        solver->add(ac == bc);
    }

    // fp.close();
    gen_node_formula_z3(a, windowDepth, solver, ctx);
    gen_node_formula_z3(b, windowDepth, solver, ctx);

    // fp.open("valid.smt2", std::ios::app);
    // fp << "(check-sat)\n";
    // fp.close();

    // std::string command = "boolector valid.smt2 > valid_output";
    // system(command.c_str());

    // std::ifstream FILE("valid_output");
    // std::string line;
    // getline(FILE, line);

    // if (flag == 1 || flag == -1)
    //     line = "unsat";

    int result = 0;
    if (solver->check() == z3::unsat)
    {
        // std::cout << solver->check() << "\n";
        if (verbose >= 3)
        {
            std::cout << a_name << " and " << b_name << " are (anti-)equivlent\n";
        }
        if (equivClass[a_name].first < equivClass[b_name].first)
        {
            if (flag == 1)
            {
                equivClass[b_name] = pair_equiv(equivClass[a_name].first, equivClass[a_name].second);
            }
            else if (flag == -1)
            {
                equivClass[b_name] = pair_equiv(equivClass[a_name].first, -equivClass[a_name].second);
            }
        }
        else
        {
            if (flag == 1)
            {
                equivClass[a_name] = pair_equiv(equivClass[b_name].first, equivClass[b_name].second);
            }
            else if (flag == -1)
            {
                equivClass[a_name] = pair_equiv(equivClass[b_name].first, -equivClass[b_name].second);
            }
        }

        result = 1;
    }
    else if (solver->check() != z3::sat)
    {
        msg("Unexpect SMT solver output!");
    }

    // free(ctx);
    delete (solver);

    return result;
}

int valid_equiv(Gate *a, Gate *b, int windowDepth, int flag, mpz_t l, mpz_t u)
{
    std::ofstream fp;
    fp.open("valid.smt2");
    fp << "(set-logic QF_BV)\n";
    node_name.clear();
    int width = NN / 3 + 1;

    char *lc = mpz_get_str(NULL, 10, l);
    char *uc = mpz_get_str(NULL, 10, u);

    if (divider_order == 2)
    {
        for (unsigned i = 0; i < NN; i++)
        {
            std::string input_name = gates[i]->get_var_name();
            node_name.insert(input_name);
        }
    }

    gen_node_declare(a, windowDepth);
    gen_node_declare(b, windowDepth);

    std::string a_name = a->get_var_name();
    std::string b_name = b->get_var_name();
    node_name.insert(a_name);
    node_name.insert(b_name);
    for (auto n : node_name)
    {
        fp << "(declare-fun " << n << "() (_ BitVec 1))\n";
    }

    /*
     * (set-logic ALL)
     * (declare-fun a0 () (_ BitVec 1))
     * (declare-fun a1 () (_ BitVec 1))
     * (declare-fun b0 () (_ BitVec 1))
     * (declare-fun b1 () (_ BitVec 1))
     * (declare-fun b2 () (_ BitVec 1))
     * (declare-fun Rn () (_ BitVec 3))
     * (declare-fun D () (_ BitVec 3))
     * (assert (= Rn (concat b2 (concat b1 b0))))
     * (assert (= D (concat (concat a1 a0) #b0)))
     * (assert (bvult Rn D))
     */
    if (divider_order == 2)
    {
        std::string concat_str_R0;
        concat_str_R0 = concat_str_divident();
        std::string concat_str_D = concat_str_divider();
        fp << "(declare-fun Rn () (_ BitVec " << 2 * width - 2 << "))\n";
        fp << "(declare-fun Dl () (_ BitVec " << 2 * width - 2 << "))\n";
        fp << "(declare-fun Du () (_ BitVec " << 2 * width - 2 << "))\n";
        fp << "(assert (= Rn " << concat_str_R0 << "))\n";
        fp << "(assert (= Dl "
           << "(bvmul " << concat_str_D << " (_ bv" << lc << " " << 2 * width - 2 << "))))\n";
        fp << "(assert (= Du "
           << "(bvmul " << concat_str_D << " (_ bv" << uc << " " << 2 * width - 2 << "))))\n";
        fp << "(assert (bvult Rn Du))\n";
        fp << "(assert (bvuge Rn Dl))\n";
    }

    fp.close();

    fp.open("valid.smt2", std::ios::app);

    if (flag == 1)
    {
        // a xor b = ((not a) and (b)) or ((a) and (not b))
        // if sat, there exists an assignment that prove a and b not equiv/anti-equiv
        // fp << "(assert (or (and (not " << a_name << ") " << b_name << ") (and " << a_name << " (not " << b_name << "))))\n";
        fp << "(assert (= #b1 (bvxor " << a_name << " " << b_name << ")))\n";
    }

    else
    {
        assert(flag == -1);
        // fp << "(assert (or (and (not " << a_name << ") "
        //    << "(not " << b_name << ")) (and " << a_name << " " << b_name << ")))\n";
        fp << "(assert (= " << a_name << " " << b_name << "))\n";
    }

    fp.close();

    gen_node_formula(a, windowDepth);
    gen_node_formula(b, windowDepth);

    fp.open("valid.smt2", std::ios::app);
    fp << "(check-sat)\n";
    fp.close();

    std::string command = "boolector valid.smt2 > valid_output";
    system(command.c_str());

    std::ifstream FILE("valid_output");
    std::string line;
    getline(FILE, line);

    // if (flag == 1 || flag == -1)
    //     line = "unsat";

    if (line == "unsat")
    {
        return 1;
    }
    else if (line != "sat")
    {
        msg("Unexpect SMT solver output!");
    }

    return 0;
}

void set_equiv(std::string aName, std::string bName, int equiv_flag)
{
    if (verbose >= 3)
    {
        {
            std::cout << aName << " and " << bName << " are (anti-)equivlent\n";
        }
    }
    if (equivClass[aName].first < equivClass[bName].first)
    {
        if (equiv_flag == 1)
        {
            equivClass[bName] = pair_equiv(equivClass[aName].first, equivClass[aName].second);
        }
        else if (equiv_flag == -1)
        {
            equivClass[bName] = pair_equiv(equivClass[aName].first, -equivClass[aName].second);
        }
    }
    else
    {
        if (equiv_flag == 1)
        {
            equivClass[aName] = pair_equiv(equivClass[bName].first, equivClass[bName].second);
        }
        else if (equiv_flag == -1)
        {
            equivClass[aName] = pair_equiv(equivClass[bName].first, -equivClass[bName].second);
        }
    }
}

/*
void valid_equiv_boolector(Gate *a, Gate *b, int windowDepth, int flag)
{
    // std::ofstream fp;
    // fp.open("valid.smt2");
    // fp << "(set-logic ALL)\n";
    Btor *btor = boolector_new();
    node_name.clear();

    gen_node_declare(a, windowDepth);
    gen_node_declare(b, windowDepth);

    std::string a_name = a->get_var_name();
    std::string b_name = b->get_var_name();
    node_name.insert(a_name);
    node_name.insert(b_name);

    BoolectorSort boolsort = boolector_bool_sort(btor);
    BoolectorNode *var;

    for (auto n : node_name)
    {
        // fp << "(declare-fun " << n << "() Bool)\n";
        var
    }
    // fp.close();
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
}
*/

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
        std::string l_name = l_gate->get_var_name();
        std::string r_name = r_gate->get_var_name();
        Gate *l_represent = gates[equivClass[l_name].first];
        Gate *r_represent = gates[equivClass[r_name].first];
        std::string l_represent_name = l_represent->get_var_name();
        std::string r_represent_name = r_represent->get_var_name();
        // node_name.insert(l_gate->get_var_name());
        // node_name.insert(r_gate->get_var_name());
        // gen_node_declare(l_gate, windowDepth - 1);
        // gen_node_declare(r_gate, windowDepth - 1);
        node_name.insert(l_represent_name);
        node_name.insert(r_represent_name);
        gen_node_declare(l_represent, windowDepth - 1);
        gen_node_declare(r_represent, windowDepth - 1);
    }
}

void gen_node_formula_z3(Gate *g, int windowDepth, z3::solver *solver, z3::context *ctx)
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
        Gate *l_represent = gates[equivClass[l_name].first];
        Gate *r_represent = gates[equivClass[r_name].first];
        std::string l_represent_name = l_represent->get_var_name();
        std::string r_represent_name = r_represent->get_var_name();

        z3::expr l_expr = name_expr.at(l_represent_name);
        z3::expr r_expr = name_expr.at(r_represent_name);

        if (equivClass[l_name].second == -1)
        { // l_represent_name = "(bvnot " + l_represent_name + ")";
            l_expr = z3::to_expr(*ctx, Z3_mk_bvnot(*ctx, l_expr));
        }
        if (equivClass[r_name].second == -1)
        {
            // r_represent_name = "(bvnot " + r_represent_name + ")";
            r_expr = z3::to_expr(*ctx, Z3_mk_bvnot(*ctx, r_expr));
        }

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

        // if (aiger_sign(l) == 0)
        //     l_literal = l_represent_name;
        // else
        //     l_literal = "(bvnot " + l_represent_name + ")";
        // if (aiger_sign(r) == 0)
        //     r_literal = r_represent_name;
        // else
        //     r_literal = "(bvnot " + r_represent_name + ")";

        if (aiger_sign(l) == 1)
            l_expr = z3::to_expr(*ctx, Z3_mk_bvnot(*ctx, l_expr));
        if (aiger_sign(r) == 1)
            r_expr = z3::to_expr(*ctx, Z3_mk_bvnot(*ctx, r_expr));
        // fprintf(f, "(assert (and %s %s))\n", l_literal.c_str(), r_literal.c_str());

        // std::ofstream fp;
        // fp.open("valid.smt2", std::ios::app);
        // fp << "(assert (= " << g_name << " (bvand " << l_literal << " " << r_literal << ")))\n";
        // fp.close();

        z3::expr current = name_expr.at(g_name);

        solver->add(current == (l_expr & r_expr));

        gen_node_formula_z3(l_represent, windowDepth - 1, solver, ctx);
        gen_node_formula_z3(r_represent, windowDepth - 1, solver, ctx);
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
        Gate *l_represent = gates[equivClass[l_name].first];
        Gate *r_represent = gates[equivClass[r_name].first];
        std::string l_represent_name = l_represent->get_var_name();
        std::string r_represent_name = r_represent->get_var_name();
        if (equivClass[l_name].second == -1)
            l_represent_name = "(bvnot " + l_represent_name + ")";
        if (equivClass[r_name].second == -1)
            r_represent_name = "(bvnot " + r_represent_name + ")";
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
            l_literal = l_represent_name;
        else
            l_literal = "(bvnot " + l_represent_name + ")";
        if (aiger_sign(r) == 0)
            r_literal = r_represent_name;
        else
            r_literal = "(bvnot " + r_represent_name + ")";
        // fprintf(f, "(assert (and %s %s))\n", l_literal.c_str(), r_literal.c_str());

        std::ofstream fp;
        fp.open("valid.smt2", std::ios::app);
        fp << "(assert (= " << g_name << " (bvand " << l_literal << " " << r_literal << ")))\n";
        fp.close();

        gen_node_formula(l_represent, windowDepth - 1);
        gen_node_formula(r_represent, windowDepth - 1);
    }
}

std::vector<Polynomial *> gen_equiv_poly_candidate(int split_num, int window, mpz_t l, mpz_t u)
{
    // print_dependency_all();

    SBIF_candidate(split_num, window, l, u);

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
            if (equiv_tag.second == 1)
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
                assert(equiv_tag.second == -1);
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

std::vector<Polynomial *> gen_equiv_poly(int split_num, int window, mpz_t l, mpz_t u)
{
    // print_dependency_all();

    SBIF(split_num, window, l, u);

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
            if (equiv_tag.second == 1)
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
                assert(equiv_tag.second == -1);
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

void print_dependency_all()
{
    std::ofstream fp;
    fp.open("dependency");
    fp.close();

    for (unsigned i = NN; i < num_gates; i++)
    {
        fp.open("dependency", std::ios::app);
        fp << "\n";
        // fp << gates[i]->get_var_name() << " ";
        fp.close();
        int depth = 0;
        max_depth = 0;
        print_dependency_gate(gates[i], depth);
    }
}

void print_dependency_gate(Gate *g, int depth)
{
    if (g->get_input())
    {
        if (depth > max_depth)
        {
            std::ofstream fp;
            fp.open("dependency", std::ios::app);
            max_depth = depth;
            fp << max_depth << " ";
            fp.close();
        }

        return;
    }
    depth++;

    std::string g_name = g->get_var_name();
    // std::cout << g->get_var_name() << " ";
    if (!g->get_output())
    {
        aiger_and *and1 = is_model_and(g->get_var_num());
        unsigned l = and1->rhs0, r = and1->rhs1;
        Gate *l_gate = gate(l), *r_gate = gate(r);
        std::string l_name = l_gate->get_var_name();
        std::string r_name = r_gate->get_var_name();
        // int l_value, r_value;
        // std::cout << aiger_sign(l) << " " << l_gate->get_var_name()
        //           << " " << aiger_sign(r) << " " << r_gate->get_var_name() << '\n';

        // l_value = (aiger_sign(l) == 1) ? 1 - gateSim[l_name].back() : gateSim[l_name].back();
        // r_value = (aiger_sign(r) == 1) ? 1 - gateSim[r_name].back() : gateSim[r_name].back();
        // gateSim[g_name].emplace_back(l_value * r_value);
        // fp << "(" << l_gate->get_var_name() << ", " << r_gate->get_var_name() << ") ";
        // fp.close();
        print_dependency_gate(l_gate, depth);
        print_dependency_gate(r_gate, depth);
    }
    else
    {
        assert(g->get_output());
        Gate *output_child_gate = g->children_front();
        std::string output_child = output_child_gate->get_var_name();
        // fp << output_child << " ";
        // fp.close();
        print_dependency_gate(output_child_gate, depth);
        // int output_child_value = (aiger_sign(slit(i - M + 1)) == 1) ? 1 - gateSim[output_child].back() : gateSim[output_child].back();
        // gateSim[g_name].emplace_back(output_child_value);
        // std::cout << aiger_sign(slit(i - M + 1)) << " " << model_output_gate->get_var_name() << "\n";
    }
}