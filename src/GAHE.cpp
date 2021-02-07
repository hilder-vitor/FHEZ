#include <vector>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>

#include "utils.h"
#include "GAHE.h"

using namespace std;
using namespace NTL;

GAHE::GAHE()  {
// Default constructor. If you use it, then you have to set everything manually!
}

GAHE::GAHE(long lambda, long gamma, long eta, long rho, long N, long log_baseG,
			ZZ t, NTL::ZZ p){

	this->lambda = lambda;
	if (N % 2 != 0){
		cout << "ERROR: N must be even (better if power of two)."
			 << "       Received N = " << N << endl;
		exit(1);
	}
	this->N = N;
	this->gamma = gamma;
	this->eta = eta;
	this->rho = rho;

	if (t < 2)
		t = ZZ(2);
	this->t = t;

	if (NumBits(p) < eta - 1)
		p = GenPrime_ZZ(eta);
	this->p = p;

	ZZ q0 = RandomBits_ZZ(gamma - eta);
	this->x0 = p * q0;

	// fmod = x^N + 1
	SetCoeff(fmod, N, 1); SetCoeff(fmod, 0, 1);

	this->sample_k_inv_k();
	// initialize b
	this->log_b = log_baseG;
	this->b = 1;
	this->b <<= log_baseG;

    double l_0 = ceil(gamma / ((double)log_baseG)); // log(2^gamma) in base b
	// Setting l sligthly bigger than l_0 because we do not perform reductions modulo x_0.
	// See Section 7 of https://eprint.iacr.org/2020/491.pdf to a discussion about 
	// the extra terms added to l_0.
    this->l = ceil(l_0 + log(N)/log(b) + log(l_0 + log(N)/log(b))/log(b)); // number of words in decomposition

	// initialize_alpha
	this->alpha = rounded_division(p, 2*t);
	
	// initialize vector g = (b^0, b^1, ..., b^(l-1))
	this->g = vector<ZZ>(this->l);
	g[0] = 1;
	for(long i = 1; i < this->l; i++){
		g[i] = g[i-1] * b;
	}

	for(long i = 0; i < this->l; i++){
		this->digits_vec.push_back(conv<ZZX>(0));
	}
}

void GAHE::sample_k_inv_k(){
	ZZ ki;
	for(long i = N-1; i >= 0; i--){
		ki = RandomBits_ZZ(eta-1);
		SetCoeff(k, i, ki);
	}

    ZZ_pPush push(p); // set the current modulus of ZZ_p to be p until the end of this function
	// inv_k = ( k^{-1} % fmod ) % p
	this->inv_k = conv<ZZX>(InvMod(  conv<ZZ_pX>(k), conv<ZZ_pX>(fmod)));
}


// Encrypts the polynomial
//       a[0] + a[1]*x + ... + a[N-1]*x^(N-1)
// into a scalar ciphertext.
NTL::ZZX GAHE::enc_scalar(const NTL::ZZX& a){
	ZZX r = sample_rx(rho, N);
	ZZX q = sample_qx(gamma, eta, N);
	ZZX msg = a; //reduce_poly_mod(msg, t); // msg = a % t
	ZZX c = p*q + r + alpha * msg;
	MulMod(c, c, k, fmod); // c = (pq + r + alpha * a) * k mod f
	reduce_poly_mod(c, x0, false); // c = c % x0 (centered = false)
	return c;
}

std::vector<NTL::ZZX> GAHE::enc_vec(const NTL::ZZX& m){
	if(deg(m) >= N){
		cout << "ERROR: trying to encrypt polynomial with degree bigger than N = " << N << endl;
		exit(1);
	}
	ZZX msg = m;
	vector<ZZX> c(l);
	for (long i = 0; i < l; i++){
		c[i] = p * sample_qx(gamma, eta, N) + sample_rx(rho, N);
		MulMod(c[i], c[i], k, fmod); // ci = (p*qi + ri) * k mod f
	}
	// at this point: c = (p * vec_q + vec_r) * k mod f
	// then, we add g * msg
	for (long i = 0; i < l; i++){
        c[i] += g[i] * msg;
	}
	// now: c = (p*vec_q + vec_r) * k + (1, b, b^2, ..., b^(l-1)) * msg   mod f
	for (long i = 0; i < l; i++)
		reduce_poly_mod(c[i], x0, false);
	// now: c = ((p*vec_q + vec_r) * k + (1, b, b^2, ..., b^(l-1)) * msg   mod f) mod x0
    return c;
}


// Receives m in [[0, 2N-1]] and encrypts x^m into a polynomial
ZZX GAHE::enc_pow_x(long m){
	ZZX msg; 
	if (m >= N){
		SetCoeff(msg, m - N, -1); // msg = -1*x^(m % N)
	}else{
		SetCoeff(msg, m, 1); // msg = x^(m % N)
	}
	return enc_scalar(msg);
}

// Receives m in [[0, 2N-1]] and encrypts x^m % (x^N+1) into a vector ciphertext
vector<ZZX> GAHE::enc_pow_x_vec(long m){
	ZZX msg; 
	if (m >= N)
		SetCoeff(msg, m - N, -1); // msg = -1*x^(m % N)
	else
		SetCoeff(msg, m, 1); // msg = x^(m % N)
    return enc_vec(msg);
}


/* Assume that c encrypts a monomial x^m, with 0 <= m <= 2N - 1,
 * then return m */
long GAHE::dec_pow_x(ZZX c){
	MulMod(c, c, inv_k, fmod);
	ZZ ai; // i-th coefficient of round((c mod p) / alpha)
	for (long i = 0; i < N; i++){
		ai = symmetric_mod(coeff(c, i), p);
		ai = rounded_division(ai, alpha);
		if (ai > 0)
			return i;
		if (ai < 0)
			return i+N;
	}
	return -1;
}

/* Assume that c encrypts a monomial x^m, then return m */
long GAHE::dec_pow_x(const vector<ZZX>& c){
	ZZX alpha_k = k * alpha;
	ZZX scalar_c = multiply(alpha_k, c);
	return dec_pow_x(scalar_c);
}


ZZX GAHE::dec(const ZZX& _c){
	ZZX c;
	MulMod(c, _c, inv_k, fmod);
	ZZX m;
	ZZ mi; // i-th coefficient of round((c*k^-1 mod p) / alpha)
	for (long i = 0; i < N; i++){
		mi = symmetric_mod(coeff(c, i), p);
		mi = rounded_division(mi, alpha);
		SetCoeff(m, i, mi);
	}
	reduce_poly_mod(m, t, false);
	return m;
}

ZZX GAHE::dec(const vector<ZZX>& c){
	ZZX alpha_k = k * alpha;
	ZZX scalar_c = multiply(alpha_k, c);
	return dec(scalar_c);
}

ZZX GAHE::multiply(const ZZX& sc,
			  const vector<ZZX>& vc){
	decompose_ZZX(digits_vec, sc, N, l, b);
	ZZX scalar_c = inner_prod(digits_vec, vc, fmod);
	return scalar_c;
}

double GAHE::get_noise(const NTL::ZZX& c, const NTL::ZZX& msg){
	ZZX u = c * inv_k;
	rem(u, u, fmod); // u = c*k^-1 % (x^N - 1) = p*q + r + alpha*msg
	u -= alpha * msg; // u = p*q + r

	double noise = 0.0;
	ZZ ri; // i-th coefficient of (u mod p)
	for (long i = 0; i < N; i++){
		ri = symmetric_mod(coeff(u, i), p);
		if (ri != 0){
			double noise_ri = log(abs(ri)) / log(2);
			if (noise_ri > noise)
				noise = noise_ri;
		}
	}
	return noise; 
}
double GAHE::get_noise(const std::vector<NTL::ZZX>& vc, const NTL::ZZX & msg){
	ZZX alpha_k = k * alpha;
	ZZX scalar_c = multiply(alpha_k, vc);
	return get_noise(scalar_c, msg);
}

// Assumes that c encrypts x^m, then compute the noise term r and
// return the logarithm of its infinity norm
double GAHE::get_noise_pow_x(const ZZX& c, long m){
	ZZX msg; 
	if (m >= N)
		SetCoeff(msg, m - N, -1); // msg = -1*x^(m % N)
	else
		SetCoeff(msg, m, 1); // msg = x^(m % N)
	return get_noise(c, msg);
}

double GAHE::get_noise_pow_x(const std::vector<NTL::ZZX>& vc, long m){
	ZZX msg; 
	if (m >= N)
		SetCoeff(msg, m - N, -1); // msg = -1*x^(m % N)
	else
		SetCoeff(msg, m, 1); // msg = x^(m % N)
	return get_noise(vc, msg);
}

std::ostream& operator<< (std::ostream &out, const GAHE& he){
	out << "GAHE: {" 
	   << "gamma: " << he.gamma
	   << ", eta: " << he.eta
	   << ", rho: " << he.rho
	   << ", N: " << he.N
	   << ", b: " << he.b
	   << ", l: " << he.l
	   << ", t: " << he.t
	   << ", p: " << he.p
	   << ", log(p): " << log(he.p)/log(2)
	   << ", log(alpha): " << log(he.alpha)/log(2)
	   << ", x0: " << he.x0
	   << "}";
	return out;
}

