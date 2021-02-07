#ifndef __GSW_LIKE_AGCD_BASED_HE__
#define __GSW_LIKE_AGCD_BASED_HE__

#include <vector>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>

// 		Randomized (Polynomial) AGCD-based Homomorphic Encryption
// 	with noise growth similar to the GSW scheme
//
//	Plaintext space is the ring Z_t[x]/<x^N + 1>.
//	We suggest to use N as a power of two so that x^N+1 is irreducible over Q.
//  When this class is used to evaluate the reresh function (bootstrapping),
// we are only interested in encrypting powers of x^2, i.e.,
// monomials like x^(2m), to represent the group Z_N.
// Therefore, to make the this usage easier, this class has functions like
// enc_pow_x, dec_pow_x, and get_noise_pow_x, which receive/return the
// exponent m instead of a polynomial.
class GAHE {

	public:
		NTL::ZZ x0;
		long N; // The polynomial ring is R = ZZ[x] / <x^N + 1>
		NTL::ZZ b; // base in which g^-1 decomposes the integers
		long log_b; // logarithm of b in base 2
		long l; // logarightm of x0 in base b
		long gamma; // bitsize of x0
		long eta; // bitsize of the prime p
		long rho; // errors sampled from ]-2^rho, 2^rho[
		long lambda; // security parameter

		NTL::ZZ alpha;  // scalar ciphertexts are of the form (p*q + r + alpha*m)*k

		NTL::ZZX fmod; // fmod = x^N + 1, defines the degree of the plaintext space
		NTL::ZZ t;  // defines the maximum value of the coefficients of plaintexts

	//private:
		NTL::ZZ p;
		NTL::ZZX k; // degree-(N-1) random polynomial with coefficients in Z_x0
		NTL::ZZX inv_k; // inverse of k in Z_p[x] / <x^N + 1>
		
		std::vector<NTL::ZZ> g; // vector of powers: (b^0, b^1, ..., b^(l-1))

		std::vector<NTL::ZZX> digits_vec; // to store g^-1 of c when multiplying c by a vector ciphertext v

	GAHE();

	GAHE(long lambda, long gamma, long eta, long rho, long N, long log_baseG,
			NTL::ZZ t = NTL::ZZ(2), NTL::ZZ p = NTL::ZZ(1));


	// sample two polynomials k and inv_k such that
	//   (k * inv_k % (x^N+1)) % p = 1
	void sample_k_inv_k();

	// Encrypts x^2m into a scalar ciphertext
	NTL::ZZX enc_pow_x(long m);
	// Encrypts x^2m into a vector ciphertext
	std::vector<NTL::ZZX> enc_pow_x_vec(long m);

	// Encrypts the polynomial
	//       a[0] + a[1]*x + ... + a[N-1]*x^(N-1)
	// into a scalar ciphertext.
	NTL::ZZX enc_scalar(const NTL::ZZX& a);
	// Encrypts the given polynomial into a vector ciphertext
	std::vector<NTL::ZZX> enc_vec(const NTL::ZZX& m);

	// return a value m in [[0, N-1]]
	long dec_pow_x(NTL::ZZX c);
	long dec_pow_x(const std::vector<NTL::ZZX>& c);

	NTL::ZZX dec(const NTL::ZZX& c);

	NTL::ZZX dec(const std::vector<NTL::ZZX>& c);

	// sc is a scalar encryption of some polynomial m1.
	// vc is a vector encryption of some polynomial m2.
	//    This function homomorphically multiplies sc by vc and returns
	// a scalar ciphertext encrypting m1*m2 mod (x^N+1)
	NTL::ZZX multiply(const NTL::ZZX& sc,
				  const std::vector<NTL::ZZX>& vc);

	// Computes the noise term r, represents it as a coefficient vector,
	// computes the infinity norm of this vector and return the log_2 of it.
	double get_noise(const NTL::ZZX& c, const NTL::ZZX& msg);
	double get_noise(const std::vector<NTL::ZZX>& vc, const NTL::ZZX& msg);
	
	// Assumes that c encrypts x^msg, then compute the noise term r and
	// return the logarithm of its infinity norm
	double get_noise_pow_x(const NTL::ZZX& c, long msg);
	double get_noise_pow_x(const std::vector<NTL::ZZX>& c, long msg);
};
	
std::ostream& operator<< (std::ostream &out, const GAHE& he);

#endif
