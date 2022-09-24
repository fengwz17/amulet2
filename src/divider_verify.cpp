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

  print_circuit_poly(f1);

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
  slicing_elim_time = process_time();

  const Polynomial *rem = reduce_divider(f2);
  if (rem && !rem->is_constant_zero_poly())
  {
    msg("INCORRECT MULTIPLIER");
    msg("");

    if (inp_f && gen_witness_divider)
    {
      msg("REMAINDER IS");
      fputs("[amulet2] ", stdout);
      rem->print(stdout);
      msg("");
      generate_witness(rem, inp_f);
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

void print_circuit_poly(FILE *file)
{
  for (unsigned i = NN; i < num_gates; i++)
  {
    Polynomial *p = gen_gate_constraint(i);
    assert(p);

    fprintf(file, "f%i = ", p->get_idx());
    p->print(file);
    delete (p);
  }
}

Polynomial *print_spec_poly_divider(FILE *file)
{
  mpz_t coeff;
  mpz_init(coeff);

  // dividend
  for (int i = NN / 2 - 1; i >= 0; i--)
  {
    const Var *a = gates[a0 + i * ainc]->get_var();
    mpz_pow_ui(coeff, base, i);
    mpz_neg(coeff, coeff);
    // if (i == static_cast<int>(NN - 1) && signed_mult)
    //   mpz_neg(coeff, coeff);

    Term *t1 = new_term(a, 0);
    Monomial *m1 = new Monomial(coeff, t1);
    push_mstack_end(m1);
  }

  // partial products
  for (int i = NN / 2 - 1; i >= 0; i--)
  {
    const Var *s = gates[i + M - 1]->get_var();

    for (int j = NN / 2 - 1; j >= 0; j--)
    {
      const Var *b = gates[b0 + j * binc]->get_var();
      mpz_pow_ui(coeff, base, i + j);
      // if (i == static_cast<int>(NN / 2 - 1) && signed_mult)
      //   mpz_neg(coeff, coeff);
      // if (j == static_cast<int>(NN / 2 - 1) && signed_mult)
      //   mpz_neg(coeff, coeff);
      add_to_vstack(s);
      add_to_vstack(b);
      Term *t1 = build_term_from_stack();
      Monomial *m1 = new Monomial(coeff, t1);
      push_mstack_end(m1);
    }
  }

  // remainder
  for (int i = NN - 1; i >= NN / 2; i--)
  {
    const Var *sr = gates[i + M - 1]->get_var();
    mpz_pow_ui(coeff, base, i - NN / 2);
    Term *tr1 = new_term(sr, 0);
    Monomial *mr1 = new Monomial(coeff, tr1);
    push_mstack_end(mr1);
  }
  mpz_clear(coeff);

  Polynomial *p = build_poly();
  p->print(file);
  return p;
  // delete (p);
}

const Polynomial *reduce_divider(FILE *file)
{
  msg("");
  msg("");
  msg("started reducing");
  Polynomial *rem = 0;

  // divider spec
  Polynomial *divider_spec = print_spec_poly_divider(file);

  for (unsigned i = NN; i < num_gates; i++)
  {
    Polynomial *p = gen_gate_constraint(i);
    // assert(p);

    // fprintf(file, "f%i = ", p->get_idx());
    // p->print(file);
    // delete (p);
    rem
  }

  return rem;
}