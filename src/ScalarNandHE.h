
#include <NTL/ZZ.h>

/**  ACD-based homomorphic scheme over {0, 1}.
*  We can perform one one OR, or one AND, or one NAND, or several NOT gates,
* then, we are supposed to bootstrap
**/
class ScalarNandHE{
	public:
		long gamma; // bitsize of ciphertexts
		long eta; // bitsize of the prime p
		long rho; // errors sampled from ]]-2^rho, 2^rho[[
		long mu; // first mu bits of ciphertexts are set to zero

		NTL::ZZ two_to_mu; // 2^mu

		NTL::ZZ alpha;  // shift applied to message before encrypting 

		NTL::ZZ p5over8; // encryption of the form pq + r + 5*p/8 used in NAND gate

	//private:
		NTL::ZZ p;


	ScalarNandHE(long gamma, long eta, long rho, long mu = 0, NTL::ZZ p = NTL::ZZ(1));


    /** Level-1 ciphertexts: c = pq + r + alpha*m with |r| < alpha/8 and m in {0,1}
    *   Level-2 ciphertexts: c = pq + r + 2*alpha*m with |r| < 3*alpha/4 and m in {0,1}
	**/
	NTL::ZZ enc(long m, long level=1);

    long dec(NTL::ZZ c, long level = 1);

    double get_noise(NTL::ZZ c, long msg, long level);

    // Homomorphic operations
	NTL::ZZ hom_nand(NTL::ZZ c1, NTL::ZZ c2); // level-1 c1 and c2 |--> level-2 c
};

// defining a way to print this class
std::ostream& operator<< (std::ostream &out, const ScalarNandHE& he);

