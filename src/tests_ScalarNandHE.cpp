#include "ScalarNandHE.h"
#include "utils.h"

bool check_mu_bits_zero(ZZ c, long mu){
	// check if first mu bits are zero	
	for (long j = 0; j < mu; j++){
		if ( (c % 2) != 0) {
			cout << "ERROR:" << endl;
			cout << j << "-th bit of c is not zero" << endl;
			cout << "mu = " << mu << endl;
			return false;
		}
		c /= 2;
	}
	return true;
}


bool test_enc_dec(ScalarNandHE& he, long NUMBER_TESTS = 75){
	for (long level = 1; level <= 2; level++){
		for (long i = 0; i < NUMBER_TESTS; i++){
    	    long m = rand_bit();
        	ZZ c = he.enc(m, level);
			if (NumBits(c) > he.gamma){
				cout << "ERROR:" << endl;
				cout << "NumBits(c) = " << NumBits(c) 
					 << " > gamma = "<<he.gamma<< endl;
				return false;
			}

	        if ((he.dec(c, level) != m) || (he.get_noise(c, m, level) > he.rho)){
				cout << "ERROR:" << endl;
				cout << "m = " << m << "    | "
					 << "he.dec(c, level="<<level<<") = " << he.dec(c, level) << endl;
				cout << "rho = " << he.rho << " | "
					 << "noise(c) = " << he.get_noise(c, m, level) << endl;
				return false;
			}

			if (false == check_mu_bits_zero(c, he.mu))
				return false;
		}
	}
    return true;
}

bool test_NAND(ScalarNandHE& he, long NUMBER_TESTS = 50){
	ZZ bound_noise = ZZ(1);
	bound_noise <<= (he.rho + 1);
	bound_noise += he.alpha/2; // bound_noise = 2^rho + alpha / 2
	double log_bound_noise = log(bound_noise) / log(2);
    for(long i = 0; i < NUMBER_TESTS; i++){
        long m1 = rand_bit();
        long m2 = rand_bit();
        ZZ c1 = he.enc(m1);
        ZZ c2 = he.enc(m2);
        ZZ cnand = he.hom_nand(c1, c2);
		long dec_nand = he.dec(cnand, 2);
		double noise_nand = he.get_noise(cnand, 1-m1*m2, 2);
 
		if ((dec_nand != (1-m1*m2)) || (noise_nand > log_bound_noise)){
				cout << "ERROR:" << endl;
				cout << "m1 = " << m1 << endl;
				cout << "m2 = " << m2 << endl;
				cout << "NAND(m1*m2) = " << 1-m1*m2 << endl;
				cout << "he.dec(cnand, level="<<2<<") = " << dec_nand << endl;
				cout << "rho = " << he.rho << endl;
				cout << "noise(cnand) = " << noise_nand << endl;
				cout << "bound noise = " << log_bound_noise << endl;
				return false;
		}
		if (false == check_mu_bits_zero(cnand, he.mu))
				return false;

		if (log(cnand)/log(2) >= he.gamma + 1){
			cout << "ERROR:" << endl;
			cout << "bit length of cnand = "
				 << log(cnand)/log(2)
				 << " > gamma + 1 = "
				 << he.gamma + 1 << endl;
			return false;			
		}
	}

    return true;
}

int main(){
	long sec_level = 100;
    long eta = sec_level + 5;
    long rho = sec_level;
    long gamma = ceil( sec_level * (eta - rho) * (eta - rho) / log(sec_level));
	long mu = rho;

    ScalarNandHE he = ScalarNandHE(gamma, eta, rho, mu);

    cout << he << endl;

    if (test_enc_dec(he)){
        cout << "Encyption and Decryption: OK" << endl;
	}else{
		cout << "ERROR:   Encyption and Decryption" << endl;
		return 1;
	}

	if (test_NAND(he)){
        cout << "Homomorphic NAND OK" << endl;
	}else{
        cout << "ERROR: Homomorphic NAND" << endl;
		return 1;
	}

	return 0;
}
