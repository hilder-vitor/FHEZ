
#include "NTT.h"
#include "utils.h"
#include <iostream>

using namespace NTL;
using namespace std;


NTT::NTT(long int log_dimension, NTL::ZZ modulus, NTL::ZZ root){
	this->logN = log_dimension;
	this->N = 1 << logN;
	this->p = modulus;
	this->omega = root;
	this->inv_omega = InvMod(omega, p); 

	this->powers_omega.SetLength(N/2);
	this->powers_inv_omega.SetLength(N/2);
	powers_omega[0] = ZZ(1);
	powers_inv_omega[0] = ZZ(1);

	// precompute powers of omega and omega^-1
	for (long int i = 1; i < N / 2; i++){
		powers_omega[i] = (powers_omega[i-1] * omega) % p;
		powers_inv_omega[i] = (powers_inv_omega[i-1] * inv_omega) % p;
	}

	this->invN = InvMod(ZZ(N), p);

	this->bit_rev_perm = vector<long int>(N);

	// precompute bit-reverse permutation
	for (long int i = 0; i < N; i++) {
		bit_rev_perm[i] = reverse_bits(i, logN);
	}

}

void NTT::transform(vec_ZZ &vec) {
	if (vec.length() != this->N){
		cout << "ERROR: vec.length() = " << vec.length()
			 << " but it must be equal to N = 2**" << logN << endl;
		exit(1);
	}
	apply_bit_reverse_permutation(vec);
	low_level_transform(vec, powers_omega);
}

void NTT::inv_transform(vec_ZZ &vec) {
	if (vec.length() != this->N){
		cout << "ERROR: vec.length() = " << vec.length()
			 << " but it must be equal to N = 2**" << logN << endl;
		exit(1);
	}
	apply_bit_reverse_permutation(vec);
	low_level_transform(vec, powers_inv_omega);

	for (long int i = 0; i < N; i++){
		vec[i] *= invN;
		vec[i] %= p;
	}
}

void NTT::apply_bit_reverse_permutation(NTL::vec_ZZ& u){
	for (long int i = 0; i < N; i++) {
		long int j = this->bit_rev_perm[i];
		if (j > i){
			ZZ tmp = u[i];
			u[i] = u[j];
			u[j] = tmp;
		}
	}
}

// Traditional radix-2 transform implemented in-place and over Z/pZ
void NTT::low_level_transform(vec_ZZ &vec, const vec_ZZ& powers) {
	
	ZZ left, right;
	for (long int size = 2; size <= N; size *= 2) { // log(N) steps
		long int halfsize = size / 2;
		long int tablestep = N / size;
		for (long int i = 0; i < N; i += size) {
			for (long int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {

				int l = j + halfsize;
				left = vec[j];
				right = vec[j + halfsize] * powers[k];
					
				vec[j] = (left + right) % p;
				vec[l] = (left - right) % p;
			}
		}
	}
}

