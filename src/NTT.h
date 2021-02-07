
#ifndef __NTT__
#define __NTT__

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

#include <vector>


/**
 *	Performs discret Fourier transform of dimension N (power of two) over Z/pZ (for prime p).
 *
 *	Notice that this p is not the same as the secret p defined in the AGCD problem.
 **/
class NTT {

	public:

		long int N;    // dimension
		long int logN; // log of N in base 2
		NTL::ZZ p;     // modulus. Must be a prime of the form N*k + 1. 
		NTL::ZZ omega; // primitive N-th root of unity on Z/pZ (order(omega) = N in Z/pZ)
		NTL::ZZ inv_omega; // omega^-1 in Z/pZ
		NTL::vec_ZZ powers_omega;  // precomputed omega^i % q for 0 <= i < N/2
		NTL::vec_ZZ powers_inv_omega; // precomputed omega^-i % q for 0 <= i < N/2

		NTL::ZZ invN; // N^-1 modulo p

		std::vector<long int> bit_rev_perm; // precomputed bit-reserve permutation


	/** 
	 *	   The dimension N must be a power of two, therefore, the constructor
	 *	receives the logarithm of the dimension in base two and set
	 *	N = 2**log_dimension.
	 */
	NTT(long int log_dimension, NTL::ZZ modulus, NTL::ZZ root);

	/** 
	 * 		Receives an N-dimensional u and compute its (forward) NTT (number 
	 * theoretic transform, that is, the discret Fourier transform over Z/pZ).
	 *		The computation is done in place, thus, the result is stored in u.
	 */
	void transform(NTL::vec_ZZ& u);
	
	/**
	 *		This is the inverse of the function transform described above, that
	 *	is, this function performs (inplace) the backwards NTT transform.
	 **/
	void inv_transform(NTL::vec_ZZ &vec);


	/** ---------------------------------------------
	 * 		Low-level functions for NTT  
	 * ---------------------------------------------- **/

	void apply_bit_reverse_permutation(NTL::vec_ZZ& u);

	/**
	 *		Let omega be the N-th primitive root of unity in Z/pZ. When this 
	 *	function is used to perform a forward transformation, one must use
	 *			powers = (1, omega, omega^2, ..., omega^(N/2 - 1) )
	 *		For the inverse transformation, one must use
	 *			powers = (1, omega^-1, omega^-2, ..., omega^-(N/2 - 1) )
	 */
	void low_level_transform(NTL::vec_ZZ &vec, const NTL::vec_ZZ& powers);
	
};

#endif
