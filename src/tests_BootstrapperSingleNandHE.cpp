#include "utils.h"
#include "BootstrapperSingleNandHE.h"

using namespace std;
using namespace NTL;

using namespace std::chrono; // to measure execution times


bool test_bootstrap_SingleNandHE(long lambda){

	double root_hermite_factor;
	if (100 == lambda)
		root_hermite_factor = 1.0064;
	else
		// XXX for other security levels, set a suitable root Hermite factor
		root_hermite_factor = 1.04;


	//################################################################################
	//### Parameters of base scheme (the one that will be bootstrapped)
	long _k_ = 5;
	long rho1 = lambda;
	long eta1 = lambda + _k_;
	long gamma1 = ceil( (eta1 - rho1) * (eta1 - rho1) / (4 * log(root_hermite_factor)/log(2)));
	if (gamma1 < 2*eta1)
		gamma1 = 2 * eta1;
	long mu = rho1 - 5;

	ScalarNandHE base_he = ScalarNandHE(gamma1, eta1, rho1, mu);

	// -----------------------------------------------------------------------
	// Parameters defining how the ciphertexts of the base scheme are decomposed 
	// B is the decomposition base.
	//  Bigger B implies larger bootstrapping keys but faster refreshing 
	long B = 1 << 6; // Decomposition base
	long L = ceil((gamma1 + 1) * log(2)/log(B)); // number of words

	// ------------------------------------------------------------------------
	// Parameters of the GSW scheme
	
	double Delta = floor((L - lambda * log(2)/log(B)) / 6);
	long N = ceil(16 * Delta);
	long log_N = ceil(log(N)/log(2));
	N = 1 << log_N;
	cout << "Delta = " << Delta << endl;

	long eta2 = lambda;
	long rho2 = floor(eta2 - sqrt(8*eta2*N*log(root_hermite_factor) / log(2)));
	cout << "rho2 = " << rho2 << endl;
	long gamma2 = ceil( (eta2 - rho2) * (eta2 - rho2) / (4 * N * log(root_hermite_factor)/log(2)));
	if (gamma2 < 2*eta2)
		gamma2 = 2 * eta2;

	long log_b = rho1 - rho2 - eta1 + eta2 - log(N * N * N * L * gamma2)/2 - 1;
	long b = 1 << log_b; // 2^(log_b) is the base in which GSW ciphertexts are decomposed during multiplication

	// instantiate GSW-like polynomial scheme
	GAHE gsw = GAHE(lambda, gamma2, eta2, rho2, N, log_b, ZZ(4));
	cout << "log(gsw.b) = " << log_b << endl;

	// ----------------------------------------------------------------------
	// Generating bootstrapping keys
	BootstrapperSingleNandHE bootstrapper = BootstrapperSingleNandHE(base_he, gsw, Delta, B);

	cout << endl << bootstrapper << endl << endl;
	cout << "Total size of refreshing keys: " << bootstrapper.get_size_bootstrapping_key_in_MB() << "MB" << endl;
	cout << "Total size of truncated refreshing keys: " << bootstrapper.get_size_trunc_bootstrapping_key_in_MB() << "MB" << endl;

	cout << "Generating bootstrapping keys" << endl;
	bootstrapper.genBootstrappingKeys();


	// ----------------------------------------------------------------------
	// Testing refresh
	long NUMBER_TESTS = 20;
	double avg_refresh_time = 0;
	high_resolution_clock::time_point t_before;
    high_resolution_clock::time_point t_after;

	for(long attempt = 1; attempt <= NUMBER_TESTS; attempt++){
	    cout << endl;
    	cout << "attempt = " << attempt << endl;
		long m1 = rand_bit();
		cout << "m1 = " << m1 << endl;
		long m2 = rand_bit();
		long msg = 1 - m1 * m2;
		cout << "m2 = " << m2 << endl;

		ZZ c1 = base_he.enc(m1);
		ZZ c2 = base_he.enc(m2);
		ZZ c_nand = base_he.hom_nand(c1, c2);

		// decrypt and test noise to debug
		long dec_c_nand = base_he.dec(c_nand, 2);
		if (dec_c_nand != msg){
			cout << "ERROR:" << endl;
			cout << dec_c_nand << " = dec_c_nand != msg = " << msg << endl;
			return false;
		}
		double noise_c_nand = base_he.get_noise(c_nand, msg, 1);

	// ----------------------------------------------------------------------
    // Start to refresh the ciphertext
		t_before = high_resolution_clock::now();

   		ZZ refreshed = bootstrapper.refresh(c_nand);
		
		t_after = high_resolution_clock::now();
    // Finish refreshing
		double t_duration = duration_cast<milliseconds>(t_after - t_before).count();
		cout << t_duration << " milliseconds to refresh a ciphertext." << endl;
		avg_refresh_time += t_duration / NUMBER_TESTS;

		// --------------------------------------------------------
		// Test if refreshed ciphertext is correct
	    long dec_refreshed = base_he.dec(refreshed, 1); // level one
		if (dec_refreshed != msg){
			cout << "ERROR:" << endl;
			cout << "dec_refreshed = " << dec_refreshed << endl;
			cout << "original msg = " << msg << endl;
			return false;
		}else{
		    cout << "Noise before refreshing: " << (base_he.get_noise(c_nand, msg, 2)) << endl;
    		cout << "Noise of refreshed ciphertext: " << (base_he.get_noise(refreshed, msg, 1)) << endl;
		    cout << "...log(refreshed, 2) = " << (log(abs(refreshed))/log(2)) << endl;
		}

    }
	cout << endl;
	cout << "AVG refreshing time: "<< avg_refresh_time << " milliseconds." << endl; 
	cout << endl;
	return true;
}


int main(){
	long lambda = 100; // Security Level
	if (test_bootstrap_SingleNandHE(lambda)){
		cout << "test_bootstrap_SingleNandHE(" << lambda<< "):  OK" << endl;
		return 0;
	}
	cout << "ERROR test_bootstrap_SingleNandHE(" << lambda<< ")" << endl;
	return 1;
}
