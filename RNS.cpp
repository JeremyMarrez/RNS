#include "RNS.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/LLL.h>
#include <assert.h>

/* precomputations */
std::array<value_t, h2> b2i_over_B2;
std::array<value_t, h1> m_1;
std::array<value_t, h2> B2_over_b2i_B1;
std::array<std::array<value_t, h1>, h2> m_2;
std::array<value_t, h2> b2i_inv_sk;
value_t B1_inv_sk;
value_t B2_inv_sk;
std::array<value_t, h1> m_sk;
std::array<std::array<value_t, h2>, h1> B2_over_b2i_mod_b1i;
std::array<value_t, h1> B2_mod_b1i;

/* channel arithmetic */
std::array<value_t, h1> c1;
std::array<value_t, h2> c2;
value_t mask;
value_t mask_sk;
  
value_t reduce (greater_value_t x, value_t c, value_t p);
value_t addmod (greater_value_t a, value_t b, value_t p);
value_t submod (value_t a, value_t b, value_t p);
value_t mulmod (greater_value_t a, value_t b, value_t c, value_t p);

template<size_t h>
std::array<value_t, h>
rns_pol (const NTL::ZZ &x,
	 const std::array<value_t, h> &ps)
{
  std::array<value_t, h> x_red;

  for (size_t i = 0; i < h; i++) {
    NTL::ZZ pi (ps[i]);
    NTL::ZZ x_red_i = x % pi;
    x_red[i] = NTL::conv<value_t> (x_red_i);
  }

  return x_red;
}

template<size_t h>
NTL::ZZ
rns_pol_inv (const std::array<value_t, h> &x1,
	     const std::array<value_t, h> &ps)
{
  NTL::ZZ X;
  NTL::ZZ P;

  if (x1[0] > (ps[0]/2)) {
      X = NTL::ZZ (x1[0]) - NTL::ZZ (ps[0]);
  } else {
    X = NTL::ZZ (x1[0]);
  }
  P = ps[0];

  for (size_t i = 1; i < h; i++) {
    NTL::CRT (X, P, NTL::ZZ (x1[i]), NTL::ZZ (ps[i]));
  }

  return X;
}

RNS::RNS ()
{
  /* channel arithmetic */
  for (size_t i = 0; i < h1; i++) {
    c1[i] = (((greater_value_t)1)<<logw) - b1[i];
  }
  for (size_t i = 0; i < h2; i++) {
    c2[i] = (((greater_value_t)1)<<logw) - b2[i];
  }

  mask = (((greater_value_t)1)<<logw) - 1;
  mask_sk = (((greater_value_t)1)<<log_bsk) - 1;

  /* precomputations */
  NTL::ZZ psk (1);
  psk <<= log_bsk;

  NTL::ZZ B2 (1);
  for (size_t i = 0; i < h2; i++) {
    B2 *= NTL::ZZ (b2[i]);
  }
  NTL::ZZ B1 (1);
  for (size_t i = 0; i < h1; i++) {
    B1 *= NTL::ZZ (b1[i]);
  }
  
  for (size_t i = 0; i < h2; i++) {
    NTL::ZZ B2_over_b2i (B2 / NTL::ZZ (b2[i]));
    b2i_over_B2[i] = NTL::conv<value_t>
      (NTL::InvMod (B2_over_b2i % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));
    value_t B1inv = NTL::conv<value_t>
      (NTL::InvMod (B1 % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));
    B2_over_b2i_B1[i] = NTL::conv<value_t>
      ((B2_over_b2i * NTL::ZZ (B1inv)) % NTL::ZZ (b2[i]));

    for (size_t j = 0; j < h1; j++) {
      m_2[i][j] = NTL::conv<value_t>
	((NTL::ZZ (b2i_over_B2[i]) * NTL::ZZ (rns_m_2[i]) *
	  NTL::InvMod (NTL::ZZ (b1[j]) % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])))
	 % NTL::ZZ (b2[i]));
    }

    for (size_t j = 0; j < h1; j++) {
      B2_over_b2i_mod_b1i[j][i] = NTL::conv<value_t>
	(B2_over_b2i % NTL::ZZ (b1[j]));
    }
  }
  
  for (size_t i = 0; i < h1; i++) {
    NTL::ZZ B1_over_b1i (B1 / NTL::ZZ (b1[i]));
    value_t b1i_over_B1 = NTL::conv<value_t>
      (NTL::InvMod (B1_over_b1i % NTL::ZZ (b1[i]), NTL::ZZ (b1[i])));
    
    m_1[i] = mulmod (rns_minv_m_1[i], b1i_over_B1, c1[i], b1[i]);

    B2_mod_b1i[i] = NTL::conv<value_t>
      (B2 % NTL::ZZ (b1[i]));
  }

  for (size_t i = 0; i < h2; i++) {
    b2i_inv_sk[i] = NTL::conv<value_t>
      (NTL::InvMod (NTL::ZZ (b2[i]) % psk, psk));
  }
  B1_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B1 % psk, psk));
  B2_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B2 % psk, psk));
  for (size_t i = 0; i < h1; i++) {
    m_sk[i] = NTL::conv<value_t>
      ((NTL::ZZ (rns_m_sk) *
	NTL::InvMod (NTL::ZZ (b1[i]) % psk, psk)) % psk);
  }
}

value_t
reduce (greater_value_t x, value_t c, value_t p)
{
  const greater_value_t bound = ((greater_value_t)1)<<logw;
  
  while (x >= bound) {
    greater_value_t x0 = x & mask;
    greater_value_t x1 = x >> logw;
    x = x0 + c * x1;
  }

  if (x >= p) x -= p;

  return x;
}

value_t
addmod (greater_value_t a, value_t b, value_t p)
{
  greater_value_t c = a + b;
  if (c >= p) c -= p;
  return c;
}

value_t
submod (value_t a, value_t b, value_t p)
{
  return (b == 0 ? a : addmod (a, p-b, p));
}

value_t
mulmod (greater_value_t a, value_t b, value_t c, value_t p)
{
  greater_value_t d = a*b;
  return reduce (d, c, p);
}

RNS_ZZq
RNS::to_ZZq (const NTL::ZZ &a)
{
  auto a_1 = rns_pol (a, b1);
  auto a_2 = rns_pol (a, b2);

  for (size_t i = 0; i < h2; i++) {
    a_2[i] = mulmod (a_2[i], b2i_over_B2[i], c2[i], b2[i]);
  }

  value_t a_sk;
  NTL::ZZ psk (1);
  psk <<= log_bsk;
  a_sk = NTL::conv<value_t> (a % psk);

  return RNS_ZZq
    {
      a_1, a_2, a_sk
    };
}

NTL::ZZ
RNS::to_ZZ (const RNS_ZZq &A)
{
  return rns_pol_inv (A.a_1, b1) % P;
}

template<size_t h>
void rns_pol_mul
(std::array<value_t, h> &zs,
 const std::array<value_t, h> &xs,
 const std::array<value_t, h> &ys,
 const std::array<value_t, h> &cs,
 const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = mulmod (xs[i], ys[i], cs[i], ps[i]);
  }
}

void rns_pol_mul_sk
(value_t &zs,
 const value_t &xs,
 const value_t &ys)
{
  zs = xs * ys;
  zs &= mask_sk;
}

template<size_t h1, size_t h2>
void
rns_pol_basis_extension (std::array<value_t, h2> &zs,
			 const std::array<value_t, h1> &xs,
			 const std::array<std::array<value_t, h1>, h2> &bs,
			 const std::array<value_t, h2> &cs2,
			 const std::array<value_t, h2> &ps2,
			 const std::array<value_t, h1> &ps1)
{
  for (size_t i = 0; i < h2; i++) {
    zs[i] = 0;
      
    for (size_t k = 0; k < h1; k++) {
      value_t pdiff = (ps1[k] > ps2[i] ? ps2[i]<<1 : ps2[i]) - ps1[k];
      value_t xskj = xs[k];
      bool negative = false;
      if (xskj > ps1[k]/2) {
	negative = true;
      }
      if (xskj > ps2[i]) {
	xskj -= ps2[i];
      }
      if (negative) {
	xskj = addmod (xskj, pdiff, ps2[i]);
      }
      value_t prod = mulmod (xskj, bs[i][k], cs2[i], ps2[i]);
      zs[i] = addmod (zs[i], prod, ps2[i]);
    }
  }
}
			 
template<size_t h1>
void
rns_pol_basis_extension_sk (value_t &zs,
			    const std::array<value_t, h1> &xs,
			    const std::array<value_t, h1> &bs,
			    const std::array<value_t, h1> &ps1)
{
  zs = 0;
      
  for (size_t k = 0; k < h1; k++) {
    value_t xskj = xs[k];
    if (xskj > ps1[k]/2) {
      xskj -= ps1[k];
    }
    value_t prod = xskj * bs[k];
    zs += prod;
  }

  zs &= mask_sk;
}

template<size_t h>
void
rns_pol_add (std::array<value_t, h> &zs,
	     const std::array<value_t, h> &xs,
	     const std::array<value_t, h> &ys,
	     const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = addmod (xs[i], ys[i], ps[i]);
  }
}

void
rns_pol_add_sk (value_t &zs,
		const value_t &xs,
		const value_t &ys)
{
  zs = xs + ys;
  zs &= mask_sk;
}

void
rns_pol_sub_sk (value_t &zs,
		const value_t &xs,
		const value_t &ys)
{
  zs = xs - ys;
  zs &= mask_sk;
}

void
rns_pol_correct (std::array<value_t, h1> &a_1,
		 const value_t &alpha,
		 const std::array<value_t, h1> B2_mod_b1i)
{
  for (size_t i  = 0; i < h1; i++) {
    value_t alpha_j = alpha;
    bool negative = false;
    if (alpha_j > (((value_t)1)<<(log_bsk-1))) {
      negative = true;
      alpha_j = -alpha_j;
      alpha_j &= mask_sk;
    }
    if (alpha_j > b1[i])
      alpha_j -= b1[i];
    value_t prod = mulmod (alpha_j, B2_mod_b1i[i], c1[i], b1[i]);
    if (negative) {
      prod = b1[i] - prod;
    }
    a_1[i] = submod (a_1[i], prod, b1[i]);
  }
}


void
RNS::mul (RNS_ZZq &c, const RNS_ZZq &a, const RNS_ZZq &b)
{
  std::array<value_t, h1> d_1;
  rns_pol_mul (d_1, a.a_1, b.a_1, c1, b1);

  value_t d_sk;
  rns_pol_mul_sk (d_sk, a.a_sk, b.a_sk);
  
  std::array<value_t, h1> xi_1;
  rns_pol_mul (xi_1, d_1, m_1, c1, b1);
  
  std::array<value_t, h2> q_2;
  rns_pol_basis_extension<h1, h2> (q_2, xi_1, m_2, c2, b2, b1);
  
  value_t q_sk;
  rns_pol_basis_extension_sk (q_sk, xi_1, m_sk, b1);
  
  rns_pol_mul (c.xi_2, a.xi_2, b.xi_2, c2, b2);
  rns_pol_mul (c.xi_2, c.xi_2, B2_over_b2i_B1, c2, b2);
  rns_pol_add (c.xi_2, c.xi_2, q_2, b2);

  rns_pol_mul_sk (d_sk, d_sk, B1_inv_sk);
  rns_pol_add_sk (c.a_sk, q_sk, d_sk);

  value_t r_over_B2;
  rns_pol_mul_sk (r_over_B2, c.a_sk, B2_inv_sk);
  value_t alpha_sk;
  rns_pol_basis_extension_sk (alpha_sk, c.xi_2, b2i_inv_sk, b2);
  rns_pol_sub_sk (alpha_sk, alpha_sk, r_over_B2);
  
  rns_pol_basis_extension (c.a_1, c.xi_2, B2_over_b2i_mod_b1i, c1, b1, b2);
  rns_pol_correct (c.a_1, alpha_sk, B2_mod_b1i);
}

RNS_ZZq
RNS::mul (const RNS_ZZq &a, const RNS_ZZq &b)
{
  RNS_ZZq c;
  mul (c, a, b);
  return c;
}
