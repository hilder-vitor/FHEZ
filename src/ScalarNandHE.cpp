
#include "ScalarNandHE.h"
#include "utils.h"

using namespace std;
using namespace NTL;

ScalarNandHE::ScalarNandHE(long gamma, long eta, long rho, long mu, ZZ p){
	this->gamma = gamma;
	this->eta = eta;
	this->rho = rho;
	if (mu >= rho)
		mu = rho - 1;
	this->mu = mu;
	two_to_mu = ZZ(1) << mu; // 2^mu

	if (NumBits(p) < eta - 1)
		p = GenPrime_ZZ(eta);
	this->p = p;

    this->alpha = rounded_division(p, ZZ(4));

	ZZ r, q;
	r = sample_r(rho-2); 
	do {
		q = sample_q(gamma, eta);
	}while(log(q)/log(2) < gamma - eta - 1);
	this->p5over8 = p*q + r + rounded_division(5*p, ZZ(8));
	this->p5over8 -= (p5over8 % two_to_mu);
}


/** Level-1 ciphertexts: c = pq + r + p*m/4 with |r| < p/16 and m in {0,1}
*   Level-2 ciphertexts: c = pq + r + p*m/2 with |r| < p/4 and m in {0,1}
**/
ZZ ScalarNandHE::enc(long m, long level){
   ZZ r, q, c;
   r = sample_r(rho-1); 
   q = sample_q(gamma-2, eta);
   if (1 == level)
	   c = p*q + r + alpha * m;
   else
	   c = p*q + r + 2 * alpha * m;
   return c - (c % two_to_mu);
}

long ScalarNandHE::dec(ZZ c, long level){
	if (1 == level){
		c = symmetric_mod(c, p);
	    return (conv<long>(rounded_division(4 * c, p)) % 2);
	}else{
		ZZ r = symmetric_mod(c - p / 2, p); // outputs r if c encrypts 1
		return (4*abs(r) < p);
	}
}

double ScalarNandHE::get_noise(ZZ c, long msg, long level){
	ZZ mult;
	if (1 == level){
		mult = alpha;
	}else
		mult = 2*alpha;
	ZZ r = symmetric_mod(c - mult * msg, p);
	if (0 == r)
		return 0;
	else
		return log(abs(r)) / log(2);
}

ZZ ScalarNandHE::hom_nand(ZZ c1, ZZ c2){ // level-1 c1 and c2 |--> level-2 c
	return this->p5over8 - c1 - c2;
}

std::ostream& operator<< (std::ostream &out, const ScalarNandHE& he){
	out << "Scalar Single NAND HE: {"
	   << "gamma: " << he.gamma
	   << ", eta: " << he.eta
	   << ", rho: " << he.rho
	   << ", mu: " << he.mu
	   << ", p: " << he.p
	   << ", log(p): " << log(he.p)/log(2)
	   << ", log(alpha): " << log(he.alpha)/log(2)
	   << "}";
	return out;
}

