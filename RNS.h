#include "Params.h"

struct RNS_ZZq {
  std::array<value_t, h1> a_1;
  std::array<value_t, h2> xi_2;
  value_t a_sk;
};

struct RNS {
  RNS ();
  RNS_ZZq to_ZZq (const NTL::ZZ &a);
  NTL::ZZ to_ZZ (const RNS_ZZq &A);
  void mul (RNS_ZZq &c, const RNS_ZZq &a, const RNS_ZZq &b);
  RNS_ZZq mul (const RNS_ZZq &a, const RNS_ZZq &b);
};
