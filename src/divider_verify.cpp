#include "divider_verify.h"
/*------------------------------------------------------------------------*/
// Global variable
bool gen_witness_divider = 1;

void dividerVerify(const char *inp_f, const char *out_f1, const char *out_f2)
{
  FILE *f1 = 0, *f2 = 0;
  if (!(f1 = fopen(out_f1, "w")))
    die("can not write output to '%s'", out_f1);

  if (!(f2 = fopen(out_f2, "w")))
    die("can not write output to '%s'", out_f2);

  print_circuit_poly_divider(f1);

  std::vector<Polynomial *> sub_poly_set = gen_substitute_poly();

  // init_slices();

  // mark_xor_chain_in_last_slice();
  init_time = process_time();

  // remove_internal_xor_gates(NULL);
  // if (!upper_half_xor_output())
  // {
  //   msg("slicing based on input cones");
  //   xor_chain = 1;
  //   slicing_non_xor();

  //   if (search_for_booth_pattern())
  //     eliminate_booth_pattern(NULL);
  //   decomposing(NULL);
  // }
  // else
  // {

  //   msg("slicing based on xor");
  //   remove_single_occs_gates(NULL);

  //   slicing_xor();
  //   remove_slice_minus_one_gates(NULL);
  // }
  // slicing_elim_time = process_time();

  // divider spec
  Polynomial *divider_spec = print_spec_poly_divider(f2);

  const Polynomial *rem = reduce_base(divider_spec, sub_poly_set);

  // const Polynomial *rem = reduce_divider(f2);
  if (rem && !rem->is_constant_zero_poly())
  {
    msg("INCORRECT MULTIPLIER");
    msg("");

    // if (inp_f && gen_witness_divider)
    if (inp_f)
    {
      msg("REMAINDER IS");
      fputs("[amulet2] ", stdout);
      rem->print(stdout);
      msg("");
      // generate_witness(rem, inp_f);
    }
  }
  else
  {
    msg("");
    msg("CORRECT MULTIPLIER");
    msg("");
    msg("writing gate constraints to '%s' ", out_f1);
    msg("writing specification to '%s'    ", out_f2);
  }
  delete (rem);

  reduction_time = process_time();

  fclose(f1);
  fclose(f2);
}

void print_circuit_poly_divider(FILE *file)
{
  int max_num = 0;
  for (unsigned i = 0; i < num_gates; i++)
  {
    const Var *x = gates[i]->get_var();
    fprintf(file, "%s, ", x->get_name());
  }
  fprintf(file, "\n");
  for (unsigned i = NN; i < num_gates; i++)
  {
    Polynomial *p = gen_gate_constraint(i);
    assert(p);

    fprintf(file, "poly f%i = ", p->get_idx());
    p->print(file);
    max_num = p->get_idx();
    delete (p);
  }
  fprintf(file, "ideal fI=");
  for (int i = 2; i <= max_num; i++)
  {
    fprintf(file, "f%i, ", i);
  }
}

Polynomial *print_spec_poly_divider(FILE *file)
{
  fix_order();
  mpz_t coeff;
  mpz_init(coeff);

  // remainder
  for (unsigned int i = NN - 1; i >= NN / 2; i--)
  {
    const Var *sr = gates[i + M - 1]->get_var();
    mpz_pow_ui(coeff, base, NN - 1 - i);
    Term *tr1 = new_term(sr, 0);
    Monomial *mr1 = new Monomial(coeff, tr1);
    push_mstack_end(mr1);
  }

  // partial products
  for (int i = NN / 2 - 1; i >= 0; i--)
  {
    const Var *s = gates[i + M - 1]->get_var();

    for (int j = NN / 2 - 1; j >= 0; j--)
    {
      const Var *a = gates[a0 + j * ainc]->get_var();
      mpz_pow_ui(coeff, base, NN - 2 - i - j);
      // if (i == static_cast<int>(NN / 2 - 1) && signed_mult)
      //   mpz_neg(coeff, coeff);
      // if (j == static_cast<int>(NN / 2 - 1) && signed_mult)
      //   mpz_neg(coeff, coeff);
      add_to_vstack(s);
      add_to_vstack(a);
      Term *t1 = build_term_from_stack();
      Monomial *m1 = new Monomial(coeff, t1);
      push_mstack_end(m1);
    }
  }

  // dividend
  for (int i = NN / 2 - 1; i >= 0; i--)
  {
    const Var *b = gates[b0 + i * binc]->get_var();
    mpz_pow_ui(coeff, base, NN / 2 - 1 - i);
    mpz_neg(coeff, coeff);
    // if (i == static_cast<int>(NN - 1) && signed_mult)
    //   mpz_neg(coeff, coeff);

    Term *t1 = new_term(b, 0);
    Monomial *m1 = new Monomial(coeff, t1);
    push_mstack_end(m1);
  }
  fix_order();
  mpz_clear(coeff);

  Polynomial *p = build_poly();
  p->print(file);
  return p;
  // delete (p);
}

Polynomial *divide_by_lt(const Polynomial *p, const Term *t)
{
  assert(t->size() == 1);
  const Var *v = t->get_var();
  Monomial *lm_tmp = p->get_mon(0);

  if (lm_tmp->get_term()->contains(v))
  {
    // -lt(p) / lt(pi)
    Term *t_rem = remainder(lm_tmp->get_term(), v);
    if (t_rem)
    {
      push_mstack_end(new Monomial(lm_tmp->coeff, t_rem->copy()));
    }
    else
    {
      push_mstack_end(new Monomial(lm_tmp->coeff, 0));
    }
  }
  Polynomial *pt = build_poly();
  return pt;
}

Polynomial *divide_by_lm(const Polynomial *p, const Monomial *tm)
{
  // p = 1, assume t is not constant, return 0
  if (p->is_constant_one_poly())
  {
    push_mstack_end(new Monomial(minus_one, 0));
    Polynomial *pt = build_poly();
    Polynomial *pa = add_poly(p, pt);
    delete (pt);
    return pa;
  }

  Term *t = tm->get_term();
  assert(t->size() == 1);
  Monomial *lm_tmp = p->get_mon(0);
  Term *ttmp = lm_tmp->get_term();

  while (t)
  {
    const Var *v = t->get_var();
    if (ttmp->contains(v))
    {
      // -lt(p) / lt(pi)
      // Term *t_rem = remainder(lm_tmp->get_term(), v);
      ttmp = remainder(ttmp, v);
      // if (t_rem)
    }
    else
    {
      Polynomial *pz = build_poly();
      return pz;
    }
    t = t->get_rest();
  }
  mpz_t coeff;
  mpz_init(coeff);
  mpz_mul(coeff, lm_tmp->coeff, tm->coeff);
  mpz_neg(coeff, coeff);
  if (ttmp)
  {
    // push_mstack_end(new Monomial(lm_tmp->coeff, t_rem->copy()));
    // push_mstack_end(new Monomial(lm_tmp->coeff, ttmp->copy()));
    push_mstack_end(new Monomial(coeff, ttmp->copy()));
  }
  else
  {
    // push_mstack_end(new Monomial(lm_tmp->coeff, 0));
    push_mstack_end(new Monomial(coeff, 0));
  }
  Polynomial *pt = build_poly();
  return pt;
}

const Polynomial *reduce_divider(FILE *file)
{
  msg("");
  msg("");
  msg("started reducing");

  Polynomial *rem = 0;
  Polynomial *rp = 0;
  // Polynomial *tmp;

  // divider spec
  Polynomial *divider_spec = print_spec_poly_divider(file);

  rp = divider_spec;

  while (!rp->is_constant_zero_poly())
  {
    unsigned i = NN;
    bool division = false;
    while (i < num_gates && division == false)
    {
      Gate *g = gates[i];
      if (g->get_elim())
        continue;
      if (verbose >= 4 && g->get_gate_constraint())
      {
        fputs("[amulet2] reducing by ", stdout);
        g->print_gate_constraint(stdout);
      }

      Polynomial *pi = g->get_gate_constraint();
      const Polynomial *negfactor = divide_by_lt(rp, pi->get_lt());
      if (negfactor->is_constant_zero_poly())
        i++;
      else
      {
        Polynomial *mult = multiply_poly(negfactor, pi);
        Polynomial *tmp = add_poly(rp, mult);
        delete (mult);
        division = true;
        rp = tmp;
        if (verbose >= 4)
        {
          fprintf(stdout, "[amulet2] reduced polynomial is ");
          rp->print(stdout);
        }
      }
      delete (negfactor);
    }
    if (division == false)
    {
      Monomial *lm_tmp = rp->get_mon(0);
      push_mstack_end(lm_tmp);
      Polynomial *lt_rp = build_poly();

      mpz_t negCoeff;
      mpz_init(negCoeff);
      mpz_neg(negCoeff, lm_tmp->coeff);
      push_mstack_end(new Monomial(negCoeff, lm_tmp->get_term()->copy()));
      Polynomial *neg_lt_rp = build_poly();

      // Polynomial *tmp_rem = add_poly(rem, lt_rp);
      // Polynomial *tmp_p = add_poly(rp, neg_lt_rp);
      if (rem)
      {
        rem = add_poly(rem, lt_rp);
        delete (lt_rp);
      }
      else
      {
        rem = lt_rp;
      }
      rp = add_poly(rp, neg_lt_rp);
      mpz_clear(negCoeff);
      delete (neg_lt_rp);
    }
  }

  msg("finish reducing");
  if (rem && !rem->is_constant_zero_poly())
  {
    fprintf(stdout, "[amulet2] remainder is ");
    rem->print(stdout);
  }
  else
  {
    fprintf(stdout, "[amulet2] remainder is 0;\n");
  }
  msg("");

  return rem;
}

void fix_order()
{
  // defined order needs to be fixed
  if (divider_order == 1)
  {
    Gate *tmp = gates[NN / 2];

    gates[NN / 2] = gates[NN / 2 - 1];
    gates[NN / 2 - 1] = tmp;
  }
}

const std::vector<Polynomial *> gen_substitute_poly()
{
  std::vector<Polynomial *> equiv_poly = gen_equiv_poly();
  std::vector<Polynomial *> substituted_poly;
  substituted_poly.clear();
  for (auto i : equiv_poly)
  {
    i->print(stdout);
  }
  for (unsigned i = NN; i < num_gates; i++)
  {
    Polynomial *p = gen_gate_constraint(i);
    // msg("original polynomial: ");
    // p->print(stdout);
    Polynomial *ep = reduce_base(p, equiv_poly);
    // msg("substituted polynomial: ");
    // ep->print(stdout);
    substituted_poly.emplace_back(ep);
  }
  FILE *f = 0;
  f = fopen("substituted_poly", "w");
  for (auto i : substituted_poly)
  {
    // i->print(stdout);
    i->print(f);
  }
  fclose(f);
  return substituted_poly;
}

Polynomial *reduce_base(Polynomial *p, std::vector<Polynomial *> poly_set)
{
  msg("");
  msg("");
  msg("started reducing");

  Polynomial *rem = 0;
  Polynomial *rp = 0;
  // Polynomial *tmp;

  rp = p;
  unsigned int poly_set_size = poly_set.size();

  while (!rp->is_constant_zero_poly())
  {
    unsigned i = 0;
    bool division = false;
    while (i < poly_set_size && division == false)
    {
      // Gate *g = gates[i];
      // if (g->get_elim())
      //   continue;
      Polynomial *pi = poly_set[i];
      if (verbose >= 4)
      {
        fputs("[amulet2] reducing by ", stdout);
        pi->print(stdout);
      }
      const Polynomial *negfactor = divide_by_lm(rp, pi->get_mon(0));
      if (negfactor->is_constant_zero_poly())
        i++;
      else
      {
        Polynomial *mult = multiply_poly(negfactor, pi);
        Polynomial *tmp = add_poly(rp, mult);
        delete (mult);
        division = true;
        rp = tmp;
        if (verbose >= 4)
        {
          fprintf(stdout, "[amulet2] reduced polynomial is ");
          rp->print(stdout);
        }
      }
      delete (negfactor);
    }
    if (division == false)
    {
      Monomial *lm_tmp = rp->get_mon(0);
      push_mstack_end(lm_tmp);
      Polynomial *lt_rp = build_poly();
      mpz_t negCoeff;
      mpz_init(negCoeff);
      mpz_neg(negCoeff, lm_tmp->coeff);
      Polynomial *neg_lt_rp = 0;
      if (lt_rp->is_constant_one_poly())
      {
        push_mstack_end(new Monomial(minus_one, 0));
        neg_lt_rp = build_poly();
      }
      else
      {
        push_mstack_end(new Monomial(negCoeff, lm_tmp->get_term()->copy()));
        neg_lt_rp = build_poly();
      }
      // Polynomial *tmp_rem = add_poly(rem, lt_rp);
      // Polynomial *tmp_p = add_poly(rp, neg_lt_rp);
      if (rem)
      {
        rem = add_poly(rem, lt_rp);
        delete (lt_rp);
      }
      else
      {
        rem = lt_rp;
      }
      rp = add_poly(rp, neg_lt_rp);
      mpz_clear(negCoeff);
      delete (neg_lt_rp);
    }
  }
  msg("finish reducing");
  if (rem && !rem->is_constant_zero_poly())
  {
    fprintf(stdout, "[amulet2] remainder is ");
    rem->print(stdout);
  }
  else
  {
    fprintf(stdout, "[amulet2] remainder is 0;\n");
  }
  return rem;
}

// Term *remainder_term(Term *t, Term *tv)
// {
//   Term *ttmp = t;
//   while (tv)
//   {
//     const Var *v1 = tv->get_var();
//     if (ttmp->contains(v1))
//     {
//       ttmp = remainder(ttmp, v1);
//     }
//     else
//     {
//       return t;
//     }
//     tv = tv->get_rest();
//   }
//   return ttmp;
// }
