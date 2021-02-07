#include <vector>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>

#include "GAHE.h"
#include "utils.h"
#include "PreComputedParams.hpp"

using namespace std;
using namespace NTL;

using namespace std::chrono; // to measure execution times

bool test_decompose_ZZX(GAHE& he){
	long N = he.N;
	long l = he.l;
	ZZ b = he.b;
	ZZX f;
	long NUMBER_TESTS = 10;
	vector<ZZX> words(l);
	vector<ZZX> powers(l);

	for (long i = 0; i < NUMBER_TESTS; i++){
		f = he.enc_pow_x(0); // random polynomial with gamma-bit coefficients
		decompose_ZZX(words, f, N, l, b);
		ZZX _f = ZZX(0);
		for(long j = 0; j < l; j++){
			_f += words[j] * he.g[j];	
		}
		rem(_f, _f, he.fmod);
		if (_f != f)
			return false;
	}
	return true;
}

bool test_enc_dec_pow_x_GAHE(GAHE& he){
	long NUMBER_TESTS = 30;
	for (long _i_ = 0; _i_ < NUMBER_TESTS; _i_++){
		long m = RandomBnd(2*he.N);
		ZZX c = he.enc_pow_x(m);
		
		long dec_c = he.dec_pow_x(c);
		if (dec_c != m){
			cout << "ERROR: dec_c != m" << endl;
			cout << "m = " << m << endl;
			cout << "dec_c = " << dec_c << endl;
			return false;
		}
		double noise = he.get_noise_pow_x(c, m);
		if (noise > he.rho){
			cout << "ERROR: noise(c) > rho" << endl;
			cout << "noise(c) = " << noise << endl;
			cout << "rho = " << he.rho << endl;
			cout << "m = " << m << endl;
			return false;
		}


		vector<ZZX> vc = he.enc_pow_x_vec(m); // encrypts x^2m into vc
		dec_c = he.dec_pow_x(vc);
		if (dec_c != m){
			cout << "ERROR: dec_vec(c) != m" << endl;
			cout << "m = " << m << endl;
			cout << "dec_c = " << dec_c << endl;
			return false;
		}
	}
	return true;
}

// random degree up to N-1 with coefficients between -t/2 and t/2
ZZX random_poly(long N, ZZ t){
	ZZX c;
	for (long i = N-1; i >= 0; i--){
		ZZ ci = RandomBnd(t);
		SetCoeff(c, i, ci);
	}
	return c;
}

bool test_enc_dec_GAHE(GAHE& he){
	long NUMBER_TESTS = 30;
	for (long _i_ = 0; _i_ < NUMBER_TESTS; _i_++){
		ZZX m = random_poly(he.N, he.t);
		ZZX c = he.enc_scalar(m);
		
		ZZX dec_c = he.dec(c);
		if (dec_c != m){
			cout << "ERROR: dec_c != m" << endl;
			cout << "dec_c = " << dec_c << endl;
			cout << "m = " << m << endl;
			return false;
		}
		double noise = he.get_noise(c, m);
		if (noise > he.rho){
			cout << "ERROR: noise(c) > rho" << endl;
			cout << "noise(c) = " << noise << endl;
			cout << "rho = " << he.rho << endl;
			cout << "m = " << m << endl;
			return false;
		}

		vector<ZZX> vc = he.enc_vec(m);
		dec_c = he.dec(vc);
		if (dec_c != m){
			cout << "ERROR: dec(c) != m" << endl;
			cout << "m = " << m << endl;
			cout << "dec_c = " << dec_c << endl;
			return false;
		}
	}
	return true;
}

bool test_multiplication_pow_x_GAHE(GAHE& he, long L = 4){
	long NUMBER_TESTS = 15;
	for (long _i_ = 0; _i_ < NUMBER_TESTS; _i_++){
		long mprod = _i_ % he.N;
		ZZX cprod = he.enc_pow_x(mprod);
        
		cout << "noise(c0) = " << he.get_noise_pow_x(cprod, mprod) << endl;

		string tabs = "";
		for (long j = 1; j <= L; j++){
			long mj = (_i_ + j) % he.N;
			vector<ZZX> cj = he.enc_pow_x_vec(mj);

			tabs = tabs + "   ";
            cout << tabs << "noise(c" << j << ") = "
				 << he.get_noise_pow_x(cj, mj) << endl;

            mprod = (mprod + mj) % he.N;
            cprod = he.multiply(cprod, cj);

            if (mprod != he.dec_pow_x(cprod)){
				cout << "ERROR: he.dec_scalar(cprod) = "
					 << he.dec_pow_x(cprod)
					 << " != " << mprod << endl;
				return false;
			}
			cout << tabs << "noise(cprod) = " << he.get_noise_pow_x(cprod, mprod) << endl;
			cout << tabs << "log_b(cprod) = " << log(abs(cprod[0]))/ he.log_b << endl;
		} 
		cout << endl;
	}
    return true;
}



bool test_linear_combination_Z_N_GAHE(GAHE& he, long L = 6){
	long NUMBER_TESTS = 10;

	for (long _i_ = 0; _i_ < NUMBER_TESTS; _i_++){
		vector<long> m(L);
		for (long i = 0; i < L; i++){
			m[i] = RandomBnd(he.N); 
		}
		ZZX enc_x_2sum = he.enc_pow_x(2*m[0]);
		long sum = m[0];
		cout << "noise(c0) = " << he.get_noise_pow_x(enc_x_2sum, 2*sum) << endl;
		string tabs = "";
		for (long j = 1; j < L; j++){
			vector<ZZX> cj = he.enc_pow_x_vec(2*m[j]);

			tabs = tabs + "   ";
			cout << tabs << "noise(c" << j << ") = "
				 << he.get_noise_pow_x(cj, 2*m[j]) << endl;

            enc_x_2sum = he.multiply(enc_x_2sum, cj);
			sum = (sum + m[j]) % he.N;
            long dec_sum = he.dec_pow_x(enc_x_2sum) / 2;
            if (sum != dec_sum){
				cout << "ERROR: dec_sum = "
					 << dec_sum
					 << " != " << sum << endl;
				return false;
			}
			cout << tabs << "noise(c_sum) = " << he.get_noise_pow_x(enc_x_2sum, 2*sum) << endl;
			cout << tabs << "log_b(c_sum) = " << log(abs(enc_x_2sum[0]))/ he.log_b << endl;
		}
		cout << endl;
	}
    return true;
}

bool test_mixed_multiplication_GAHE(GAHE& he, long L = 2){

	double avg_mix_prod_time = 0;
	high_resolution_clock::time_point t_before;
    high_resolution_clock::time_point t_after;

	long NUMBER_TESTS = 15;
	for (long _i_ = 0; _i_ < NUMBER_TESTS; _i_++){
		ZZX mprod = random_poly(he.N, he.t);
		ZZX cprod = he.enc_scalar(mprod);
        
		cout << "noise(c0) = " << he.get_noise(cprod, mprod) << endl;

		string tabs = "";
		for (long j = 1; j <= L; j++){
			ZZX mj = random_poly(he.N, he.t);
			vector<ZZX> cj = he.enc_vec(mj);

			tabs = tabs + "   ";
            cout << tabs << "noise(c" << j << ") = "
				 << he.get_noise(cj, mj) << endl;

			// multiply the messages in clear
			// mprod = (mprod * mj % (x^N - 1)) % t
			MulMod(mprod, mprod, mj, he.fmod);
			reduce_poly_mod(mprod, he.t, false);
			

			// ----------------------------------------------------------------------
			// Measure time to perform mixed homomorphic product
			t_before = high_resolution_clock::now();
            cprod = he.multiply(cprod, cj);
			t_after = high_resolution_clock::now();

			double t_duration = duration_cast<milliseconds>(t_after - t_before).count();
			avg_mix_prod_time += t_duration / (L * NUMBER_TESTS);


			ZZX dec_prod = he.dec(cprod);
			reduce_poly_mod(dec_prod, he.t, false);

            if (mprod != dec_prod){
				cout << "ERROR: he.dec(cprod) = "
					 << dec_prod
					 << " != " << mprod << endl;
				return false;
			}
			cout << tabs << "noise(cprod) = " << he.get_noise(cprod, mprod) << endl;
			cout << tabs << "log_b(cprod) = " << log(abs(cprod[0]))/ he.log_b << endl;
		} 
		cout << endl;
	}

	cout << endl;
	cout << "AVG time of mix_prod with N = " << he.N << ": "<< avg_mix_prod_time << " milliseconds." << endl; 
	double size_MB = he.N * he.l * he.gamma / (8.0 * 1000000);
	cout << "Size of vec ciphertext: " << size_MB << " MB" << endl;
	cout << endl;

    return true;
}



int main(){

	long lambda = 100;

	long logN = 10;

	std::map<long , PreComputedParams> map_params = params_100_private_x0;
	PreComputedParams params(map_params.at(logN));


	long eta = lambda;
	long gamma = params.gamma;
	long rho = params.rho;
	long N = params.N;
	long log_b = params.log_b;
	ZZ t = ZZ(2);

	GAHE he = GAHE(lambda, gamma, eta, rho, N, log_b, t);

	cout << he << endl;

	if (test_decompose_ZZX(he)){
		cout << "Decomposition OK" << endl;
	}else{
		cout << "ERROR: Decomposition" << endl;
		return 1;
	}
//	
//	if (test_enc_dec_pow_x_GAHE(he)){
//		cout << "Encryption and decryption of powers x^m  OK" << endl;
//	}else{
//		cout << "ERROR: Decryption or Encryption of powers x^m" << endl;
//		return 1;
//	}
//
//	if (test_enc_dec_GAHE(he)){
//		cout << "General encryption and decryption   OK" << endl;
//	}else{
//		cout << "ERROR: General decryption or encryption" << endl;
//		return 1;
//	}


//	if (test_multiplication_pow_x_GAHE(he)){
//		cout << "Multiplication scalar-vector powers of x GAHE OK" << endl;
//	}else{
//		cout << "ERROR: Multiplication scalar-vector powers of x GAHE" << endl;
//		return 1;
//	}
//
	if (test_linear_combination_Z_N_GAHE(he)){
		cout << "Linear combination modulo N using GAHE's mixed product OK" << endl;
	}else{
		cout << "ERROR: Linear combination modulo N using GAHE's mixed product" << endl;
		return 1;
	}


	if (test_mixed_multiplication_GAHE(he)){
		cout << "Multiplication scalar-vector general polynomials GAHE OK" << endl;
	}else{
		cout << "ERROR: Multiplication scalar-vector general polynomials GAHE" << endl;
		return 1;
	}

	cout << "Passed all the tests. OK" << endl;
	return 0;
}
