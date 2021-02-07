#ifndef __BOOTSTRAP_SINGLE_NAND_SCHEME__
#define __BOOTSTRAP_SINGLE_NAND_SCHEME__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

#include <vector>

#include "ScalarNandHE.h"
#include "GAHE.h"
#include "PolynomialMultiplier.h"


/**
 *   This class is used to bootstrap the "scalar" scheme based on ACD problem
 * that is implemented in the class ScalarNandHE. We call that scheme the
 * base scheme. To perform the bootstrap, the "polynomial" GSW-like scheme
 * implemented by the class GAHE is used.
 **/
class BootstrapperSingleNandHE {

	public:
		ScalarNandHE& base_scheme;
		GAHE& poly_scheme;
		// security level		
		long lambda;
		// Base in which the ciphertexts of the base scheme are decomposed during bootstrap
		long B;
		// Number of (log B)-bit words needed to decompose ciphertexts of the base scheme
		long L;
		// Bootstrapping on R := Z[x]/<x^N + 1>
		long N;
		// Upper bound to the sum of rounding errors of the hidden modulus switching
		double Delta;

		// first mu bits of ciphertexts to be refreshed are zero,
		// thus, we can ignore them. But we decompose the ciphertexts
		// in base B, hence, we can ignore mu_B := log_B(2^mu) words
		long mu_B;

		// scalar encryption of (x^2)^(shift % N)
		// In the paper, this key is named K_delta
		NTL::ZZX bk_shift;

		// bk[g][i] is a precompute NTT of a vector encryption of
		//      (x^2)^(round(g * b^i * N / p1)) % (x^N+1)
		// for 0 <= g < b, 0 <= i < l
		//
		// In the paper, these keys are named K_{g, i} while bk consists of
		// K_delta, K_{g, i}'s and ek
		std::vector< std::vector< PrecomputedVecNTT  > > bk;

		// extraction key: encrypts G * test_vector
		NTL::vec_ZZ ek;
		long gamma_ek;

		// level-1 base scheme encryption of 1/2, i.e., bar{p}*q + r + bar{p}/8
		NTL::ZZ base_scheme_p_over_8;
	
		// Auxiliar vector used to store coefficients of a polynomial at the end of refreshing
		NTL::vec_ZZ vec_coeff_refreshed;


		std::vector<NTL::ZZX> digits_vec; // to store g^-1 of c when multiplying c by a vector ciphertext v

		// multiplier is used to precompute the NTT (number theoretic transform)
		// of the bootstrapping keys to speed up the refreshing procedure
		PolynomialMultiplier* multiplier;

	BootstrapperSingleNandHE(ScalarNandHE& base_scheme, GAHE& vectorial_scheme, double Delta, long B=2);

	~BootstrapperSingleNandHE();

    void genBootstrappingKeys();

	double get_size_bootstrapping_key_in_MB();

	double get_size_trunc_bootstrapping_key_in_MB();

	NTL::ZZ refresh(NTL::ZZ c);

};
#endif

std::ostream& operator<< (std::ostream &out, const BootstrapperSingleNandHE& boot);
