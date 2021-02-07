#ifndef __MY_UTILS_FUNCS__
#define __MY_UTILS_FUNCS__

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>

#include <vector>

using namespace std;
using namespace NTL;


/* Set A = B * C mod n */
void mulMod(mat_ZZ& A, const mat_ZZ& B, const mat_ZZ& C, const ZZ& n);


/* Multiply v by the j-th colum of A and return the result reduced mod n */
ZZ mulMod_v_column(const vec_ZZ& v, const mat_ZZ& A, const ZZ& n, long j);


void reduce_vec_mod(vec_ZZ& vec, ZZ mod);

void reduce_mat_mod(mat_ZZ& mat, ZZ mod);

void reduce_poly_mod(ZZX& poly, ZZ mod, bool centered = true);

NTL::ZZ max_norm(const NTL::ZZX& poly);
NTL::ZZ max_norm(const std::vector<NTL::ZZX>& vec);

ZZ_pX inner_prod(const vector<ZZ_pX>& u, const vector<ZZ_pX>& v, const ZZ_pXModulus& f);
ZZX inner_prod(const vector<ZZX>& u, const vector<ZZX>& v, const ZZX& f);
ZZX inner_prod(const vec_ZZ& u, const vector<ZZX>& v);

bool rand_bit();

/**
 * Return  round(a/n), that is, interpret a/n as an rational number
 * then return the closest integer to it.
 */
ZZ rounded_division(const ZZ& a, const ZZ& n);

ZZX rounded_division(const ZZX& a, const ZZ& n);

/**
 * Return ceil(a/n), that is, interpret a/n as an rational number
 * then return the minimun integer bigger than a/n.
 */
ZZ ceil_division(const ZZ& a, const ZZ& n);
	

ZZ symmetric_mod(const ZZ& a, const ZZ& n); // returns a % n using the set ]-c/2, c/2] as Z/nZ

/**
 *	 Reduces the real number x modulo n, that is, returns a real value r such
 * that 0 <= r < n and a - r is an integer that is a multiple of n
 */
RR mod_on_R(const RR& x, const ZZ& n);

/**    Set words_vec to the vector (x0, x1, ..., x_{l-1})
 * representing the decomposition of x in base b.
 */
void decompose_ZZ(vec_ZZ& words_vec, const ZZ& x, long l, ZZ b);

/**
 *	Assumes that deg(f) < N and that all the coefficients of f are smaller than
 * b^l in absolute value.
 * 	Sets words_vec to the vector (y0, y1, ..., y_{l-1})
 * representing the decomposition of f in base b, that is, each yi is a
 * polynomial with degree smaller than N and coefficients in ]-b, b[, and
 * the sum y0*b^0 + y1*b^1 + ... + y_{l-1}*b^{l-1} equals f
 */
void decompose_ZZX(vector<ZZX>& words_vec, const ZZX& f, long N, long l, ZZ b);

void decompose_ZZ_pX(vector<ZZ_pX>& words_vec, const ZZ_pX& f, long N, long l, ZZ b);

/*		Receives a vector v and multiplies it by G = I.tensor_product(g),
 *	where g is the column vector (1, b, b^2, ..., b^(l-1)).
 *		The returned vector is G * v.
 *		For N = 3, for example, we have
 *			[ g 0 0 ]
 *		G = [ 0 g 0 ] \in Z^(Nl x N)
 *			[ 0 0 g ]
 *	then, notice that when we do G*v, we multiply 1, then, b, ..., then b^(l-1)
 *	by v[0], then we multiply 1, b, ..., b^(l-1), by v[1], and so on.
 */
vec_ZZ multiply_by_G(const vec_ZZ& v, long l, ZZ b);


/*
 *      Set words_vec to G^-1(vec), that is, assuming that vec is an n-dimensional vector
 * and words_vec an (n*l)-dimensional, then words_vec[0:l-1] corresponds to the base-b decomposition
 * vec[0], bits_vec[l:2*l-1] corresponds to the base-b decomposition of vec[1], and so on.
 * 		Thus l must have a value such that max(abs(vec[i])) < b^l.
 */
void invG(vec_ZZ& words_vec, const vec_ZZ& vec, long l, ZZ b);

void invG_mat(mat_ZZ& words_mat, mat_ZZ& mat, long l, ZZ b);

vec_ZZ coefficient_vector(const ZZX& f, long N);

/* Receives a polynomial fmod of degree N and a polynomial f
 * of degree smaller than N.
 * Returns an NxN matrix F such that, for any polynomial g, the following holds:
 * 		 phi(g * f % fmod) = phi(g) * F
 * where phi is a function that receives a polynomial and returns its
 * N-dimensional vector of coefficients
 **/
NTL::mat_ZZ circulant_matrix(const NTL::ZZX& f, const NTL::ZZX& fmod);

/* --------------------------------------------------------------------
* Auxiliary functions used by several encryption schemes based on AGCD
* --------------------------------------------------------------------- */

ZZ sample_r(long rho);

ZZX sample_rx(long rho, long N);

vector<ZZX> sample_vec_rx(long rho, long N, long l);

ZZ sample_q(long gamma, long eta);

ZZX sample_qx(long gamma, long eta, long N);

vector<ZZX> sample_vec_qx(long gamma, long eta, long N, long l);


std::vector<ZZX>& operator+= (std::vector<NTL::ZZX>& vec, const std::vector<NTL::ZZX>& u);
/* multiply each entry of vec by the integer t. */
std::vector<NTL::ZZX>& operator*= (std::vector<NTL::ZZX>& vec, const NTL::ZZ& t);



/**  Returns true if p is (probably) prime. 
**/
bool is_prime(const NTL::ZZ& n);


/**		Suppose that we want to multiply two polynomials f0 and f1 of degree
*	less than N. Moreover, suppose that the coefficients of f0 belong to
*	{-2^max0+1, ..., 2^max0 - 1} (and the same for f1 and max1).
*
* 		Then, this function returns a prime p = k*N + 1, for some prime k,
* 	such that, Z/pZ has a N-th primitive root of unity and we can use
* 	the NTT over Z/pZ to multiply f0 * f1 % (x^N + 1), that is, the 
* 	coefficients of the final product belong to {-p/2, ..., p/2}. */
NTL::ZZ find_modulus(int max0, int max1, long int N);

/**		Let v0 and v1 be L-dimensional vectors of polynomials of degree up to
* N-1. Moreover, considering the infinity norm, assume that	||v0|| < 2^max0 
* and ||v1|| < 2^max1, that is, the coefficients of the
* polynomials of v0 are max0-bit signed integers (and those of v1 have max1 bits).
* 		Then, this function returns a prime p = k*N + 1, for some prime k,
* such that, Z/pZ has a N-th primiitive root of unity and we can use
* the NTT over Z/pZ to compute the inner product v0 * v1 mod (x^N + 1),
* that is, the coefficients of the resulting polynomial belong to {-p/2, ..., p/2}.
*/
NTL::ZZ find_modulus(int max0, int max1, long int N, long int L);


/**
* This function assumes that N is a power of two.
*	Notice that if g^N = 1 mod p, then order(g) | N, thus, order(g) is one
* of the powers 2^1, 2^2, ..., 2^log(N) = N.
*
* 	This function returns true if g is an N-th primitive root of unity modulo p.
*/
bool is_primitive_root(NTL::ZZ g, long int N, NTL::ZZ p);

/**		Finds a N-th primitive root of unity in Z/pZ, that is, an element g
*	such that order(g) = N (the smallest k such that g^k=1 mod p is k=N).
*		This function assumes that N is a power of two and that p is a prime
* of the form N*k + 1, where k is also a prime.
* 		To find such p, you can use the function find_modulus.
* 		Thus, for all g in Z/pZ, the order of g divides N*k, which means that
* order(g) is k, or order(g) = 2^i, or order(g) = 2^i * k, for some 1 <= i <= logN.
* 		Therefore, for any g, we can test if order(g) = N or order(g^k) = N. */
NTL::ZZ find_Nth_primitive_root(long int N, NTL::ZZ p);

long reverse_bits(long int x, long int logN);


void coordinatewise_prod(const NTL::vec_ZZ u, 
						const NTL::vec_ZZ v,
						NTL::vec_ZZ& output);

/** Coordinatewise product modulo p */
void coordinatewise_prod(const NTL::vec_ZZ u, 
						const NTL::vec_ZZ v,
						NTL::vec_ZZ& output, 
						const NTL::ZZ& p);

/**
 *		Transforms u = (u_1, ..., u_n) in 
 *	 (u_1 * b^0 % p, u_2 * b^1 % p, ..., u_n * b^(n-1) % p)
 **/
void apply_powers(NTL::vec_ZZ& u, const NTL::ZZ& b, const NTL::ZZ& p);


void centralized_mod(NTL::vec_ZZ& u, const NTL::ZZ& p);

// Converts a vector of integers to a polynomial.
NTL::ZZX vec_ZZ_to_poly(const NTL::vec_ZZ& u);


/* ------------------------------------------------------
 * 			Random values generation
 * -----------------------------------------------------*/
NTL::vec_ZZ random_vec(int dimension, long int bound = 2);
	
NTL::ZZX random_poly(long int degree, long int bitlen = 2);


/**		Returns a d-dimensional vector of polynomials whose degree is up to m
 * and coefficients have absolute value smaller than 2^bitlen - 1.
 */
std::vector<NTL::ZZX> random_vec_of_polys(long int m, int d, long int bitlen = 2);

#endif
