#include "BootstrapperSingleNandHE.h"
#include "utils.h"

#include<NTL/ZZ_pX.h>

#define DEBUG false

using namespace std;
using namespace NTL;


BootstrapperSingleNandHE::BootstrapperSingleNandHE(ScalarNandHE& base_scheme, GAHE& poly_scheme, double Delta, long B)
	: base_scheme(base_scheme), poly_scheme(poly_scheme) {

	this->lambda = poly_scheme.lambda;

	this->B = B;
	this->N = poly_scheme.N;
	this->L = ceil((base_scheme.gamma + 1) / (log(B) / log(2))); // log_B(2^(gamma))
	this->N = poly_scheme.N;

	this->Delta = Delta;
	
	this->mu_B = floor(base_scheme.mu * log(2) / log(B));

	this->vec_coeff_refreshed.SetLength(N);
	
	this->poly_scheme.alpha = rounded_division(poly_scheme.p, ZZ(8));

	this->gamma_ek = base_scheme.gamma - (log(poly_scheme.N * poly_scheme.l * poly_scheme.b)/log(2));


	digits_vec = poly_scheme.digits_vec; // to store g^-1 of c when multiplying c by a vector ciphertext v


	// Initialize NTT to perform multiplications of polynomials with precomputed 
	// NTT of the bootstrapping keys
	int logN = (int)(log(N)/log(2));
	int logb = poly_scheme.log_b;
	this->multiplier = new PolynomialMultiplier(logN, poly_scheme.l, logb, poly_scheme.gamma);
}

BootstrapperSingleNandHE::~BootstrapperSingleNandHE(){
	delete multiplier;
}

void BootstrapperSingleNandHE::genBootstrappingKeys(){
	cout << "mu_B = " << mu_B << endl;

	ZZ two_to_rho = ZZ(1); two_to_rho <<= base_scheme.rho;
	long shift = conv<long>(Delta + rounded_division( (3*two_to_rho + poly_scheme.alpha/2) * N, base_scheme.p));
	cout << "shift:" << shift << endl;


	// generating the "shift key"
	this->bk_shift = poly_scheme.enc_pow_x(2*shift);
	if (poly_scheme.dec_pow_x(bk_shift)/2 != shift){
		cout << "--PROBLEM:" << endl;
		cout << "poly_scheme.dec_pow_x(bk_shift)/2 != shift" << endl;
		cout << "dec:" << poly_scheme.dec_pow_x(bk_shift)/2 << endl;
		cout << "shift:" << shift << endl;
		cout << "shift mod N:" << shift % N << endl;
		exit(1);
	}

	// generating the multiplication keys
	long exponent;
	this->bk = vector< vector< PrecomputedVecNTT > >(B);

    for(long g = 1; g < B; g++){
		if (1 == g % 10)
			cout << "g / B = " << g << " / " << B << endl;
		bk[g] = vector< PrecomputedVecNTT >(L);

		ZZ pow_B = power(ZZ(B), mu_B);
		for(long i = mu_B; i < L; i++){
			exponent = conv<long>(rounded_division( (g * pow_B * N), base_scheme.p) % N);
			bk[g][i] = multiplier->precompute_NTT(poly_scheme.enc_pow_x_vec(2*exponent));
			
			// ---- START Debug
			#if DEBUG 
				long dec_expo = poly_scheme.dec_pow_x(bk[g][i])/2;
				if ((dec_expo - exponent) % N != 0){
					cout << "PROBLEM:" << endl;
					cout << "poly_scheme.dec_pow_x(bk[" << g << "][" << i <<"])/2 != exponent" << endl;
					cout << "dec:" << endl;
					cout << dec_expo << endl;
					cout << "exponent:" << exponent << endl;
					cout << "2*exponent:" << 2*exponent << endl;
					cout << "exponent mod N:" << exponent % N << endl;
					exit(1);
				}
			#endif
			// ------ END debug
			pow_B *= B;
		}
	}


	long rho_ek = base_scheme.rho - (log(poly_scheme.N * poly_scheme.l * poly_scheme.b)/log(2));
	ZZ q0 = RandomBits_ZZ(gamma_ek - base_scheme.eta);
	ZZ x0_ek = base_scheme.p * q0; 

	vec_ZZ vec_final_interval;
	vec_final_interval.SetLength(poly_scheme.N);
	for(long i = 0; i < poly_scheme.N; i++)
		vec_final_interval[i] = 1;

	vec_final_interval = circulant_matrix(poly_scheme.inv_k, poly_scheme.fmod) * vec_final_interval;
	vec_ZZ G_invK_u = multiply_by_G(vec_final_interval, poly_scheme.l, poly_scheme.b); // G*u : N*l dimensional

	this->ek.SetLength(poly_scheme.N * poly_scheme.l);

	// in this loop, we set
	// ek = (p1*q + r + G*K^-1*u) * k in R^(N * l)
	// so that, when we multiply G^-1(c) * ek, for a c encrypting x^e,
	// where e = r + alpha * m, we get a "base scheme" encryption of
	// 		phi(x^e) * u = m
	for(long i = 0; i < ek.length(); i++){
		ek[i] = rounded_division(base_scheme.p * G_invK_u[i], poly_scheme.p);
		ek[i] += base_scheme.p*sample_q(gamma_ek, base_scheme.eta);
		ek[i] += sample_r(rho_ek);
		ek[i] = ek[i] % x0_ek;
	}

	// Set base scheme encryption of "1/2", that is,
	// c = base_scheme.p*q + r + p/8
	this->base_scheme_p_over_8 = base_scheme.p * sample_q(base_scheme.gamma - 2, base_scheme.eta)
								 + sample_r(base_scheme.rho - 10)
								 + rounded_division(base_scheme.p, ZZ(8));
}

NTL::ZZ BootstrapperSingleNandHE::refresh(NTL::ZZ c){
	vec_ZZ words;
	words.SetLength(L);
	decompose_ZZ(words, c, L, ZZ(B));

	cout << "bit length c: " << log(c) / log(2) << endl;

	ZZX refreshed = bk_shift;

	
	#if DEBUG
		ZZ pow_B = power(ZZ(B), mu_B);
		long exponent = poly_scheme.dec_pow_x(bk_shift)/2;
		cout << "shift = " << exponent << endl;
	#endif

	
    for(long i = mu_B; i < L; i++){
        long g = conv<long>(words[i]);
        if (g > 0){

			// multiply homomorphically refreshed by bk[g][i], i.e.,
			// compute v = g^-1(refreshed), then compute the inner product v * bk[g][i]
			decompose_ZZX(digits_vec, refreshed, poly_scheme.N, poly_scheme.l, poly_scheme.b);
			refreshed = multiplier->multiply(digits_vec, bk[g][i]);

			// ---- START Debug
			#if DEBUG 
				exponent += conv<long>(rounded_division( (g * pow_B * N), base_scheme.p) % N);
				exponent %= N;

				long dec_expo = poly_scheme.dec_pow_x(refreshed)/2;
				if ((dec_expo - exponent) % N != 0){
					cout << "PROBLEM:" << endl;
					cout << "i = " << i << endl;
					cout << "poly_scheme.dec_pow_x(refreshed) != exponent" << endl;
					cout << "dec:" << dec_expo << endl;
					cout << "exponent:" << exponent << endl;
					cout << "exponent mod N:" << exponent % N << endl;
					exit(1);
				}
			#endif
			// ------ END debug

		}
		#if DEBUG
		pow_B *= B;
		#endif
	}

	// ---- START Debug
	#if DEBUG 
		ZZ two_to_rho = ZZ(1) << base_scheme.rho;
		long shift = conv<long>(Delta + rounded_division( (2*two_to_rho + poly_scheme.alpha/2) * N, base_scheme.p));
		long msg = base_scheme.dec(c, 2);
		long deg_ref = poly_scheme.dec_pow_x(refreshed)/2;
		cout << "DEBUG ON" << endl;
		cout << "(x^2)^e --> e = " << deg_ref << endl;
		if (1 == msg)
        	cout << N/2 << " = N/2 <= e < N = "
				  << N << " ??  " << (N/2 <= deg_ref && deg_ref < N)
				  << endl;
	    else
			cout << "0 <= e < N/2 = " << N/2
				 << " ?? " << (deg_ref <= N/2) << endl;

	#endif
	// ------ END debug
	//
	
	cout << "log_b(norm(refreshed)) = " << log(max_norm(refreshed))/log(poly_scheme.b) << endl;
	cout << "poly_scheme.l = " << poly_scheme.l << endl;

	vec_ZZ phi_ref = coefficient_vector(refreshed, poly_scheme.N);
	words.SetLength(poly_scheme.N * poly_scheme.l);
	invG(words, phi_ref, poly_scheme.l, poly_scheme.b);

	// ---- START Debug
	#if DEBUG 
		// test if _c = p*q + r  + p/8 * (-1)^m
		ZZ _c = words * ek;
		ZZ r_alpha = symmetric_mod(_c, base_scheme.p);
		ZZ _msg = rounded_division(8 * r_alpha, base_scheme.p);
		if ((0 == msg && 1 == _msg) || (1 == msg && -1 == _msg)) {
			cout << "phi(x^(2e)) * u = (-1)^" << msg << " ...  OK" << endl;
		}else{
			cout << "ERROR:   phi(x^(2e)) * u = " << _msg << " != (-1)^m " << endl;
		}
	#endif
	// ------ END debug

	return this->base_scheme_p_over_8 - words * ek;
}

double BootstrapperSingleNandHE::get_size_bootstrapping_key_in_MB(){
	double k_shift_MB = poly_scheme.N * poly_scheme.gamma / (8.0 * 1000000);
	double each_bk_MB = poly_scheme.l * poly_scheme.N * poly_scheme.gamma / (8.0 * 1000000);
	double ek_MB = poly_scheme.N * poly_scheme.l * gamma_ek  / (8.0 * 1000000); 
	double total = k_shift_MB + (B-1) * L * each_bk_MB + ek_MB;
	return total;
}

double BootstrapperSingleNandHE::get_size_trunc_bootstrapping_key_in_MB(){
	double k_shift_MB = poly_scheme.N * poly_scheme.gamma / (8.0 * 1000000);
	double each_bk_MB = poly_scheme.l * poly_scheme.N * poly_scheme.gamma / (8.0 * 1000000);
	double ek_MB = poly_scheme.N * poly_scheme.l * gamma_ek  / (8.0 * 1000000); 
	double total = k_shift_MB + (B-1) * (L - mu_B) * each_bk_MB + ek_MB;
	return total;
}

std::ostream& operator<< (std::ostream &out, const BootstrapperSingleNandHE& boot){
	
	out << "BootstrapperSingleNandHE: {" << endl
	   << " lambda: " << boot.lambda << "," << endl
	   << " decomposition B: " << boot.B << "," << endl
	   << " decomposition L: " << boot.L << "," << endl 
	   << " base scheme: " << boot.base_scheme << "," << endl
	   << " poly scheme: " << boot.poly_scheme
	   << "}";
	return out;
}

