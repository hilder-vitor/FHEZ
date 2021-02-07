#include "utils.h"

using namespace std;
using namespace NTL;


/* Set A = B * C mod n */
void mulMod(mat_ZZ& A, const mat_ZZ& B, const mat_ZZ& C, const ZZ& n){
	mul(A, B, C);
	for (int i = 0; i < A.NumRows(); i++)
		for (int j = 0; j < A.NumCols(); j++)
			A[i][j] %= n;
}

/* Multiply v by the j-th colum of A and return the result reduced mod n */
ZZ mulMod_v_column(const vec_ZZ& v, const mat_ZZ& A, const ZZ& n, long j){
	ZZ res = ZZ(0);
	for (int i = 0; i < A.NumRows(); i++){
		res += v[i] * A[i][j];
	}
	return res % n;
}

void reduce_vec_mod(vec_ZZ& vec, ZZ mod){
	for (long i = 0; i < vec.length(); i++)
		vec[i] %= mod;
}

void reduce_mat_mod(mat_ZZ& mat, ZZ mod){
	for (long i = 0; i < mat.NumRows(); i++){
		for (long j = 0; j < mat.NumCols(); j++){
			mat[i][j] %= mod;
		}
	}
}

void reduce_poly_mod(ZZX& poly, ZZ mod, bool centered){
	ZZ ai;
	for (long i = 0; i <= deg(poly); i++){
		ai = coeff(poly, i);
		if (centered)
			SetCoeff(poly, i, symmetric_mod(ai, mod));
		else
			SetCoeff(poly, i, ai % mod);
	}
}

ZZ max_norm(const ZZX& poly){
	ZZ inf_norm = ZZ(0);
	for(long i = 0; i < deg(poly); i++){
		if (inf_norm < abs(coeff(poly, i)))
			inf_norm = abs(coeff(poly, i));
	}
	return inf_norm;
}

ZZ max_norm(const vector<ZZX>& vec){
	ZZ global_norm = ZZ(0);
	ZZ tmp_norm = ZZ(0);
	for(long i = 0; i < vec.size(); i++){
		tmp_norm = max_norm(vec[i]);
		if (global_norm < tmp_norm)
			global_norm = tmp_norm;
	}
	return global_norm;
}
// return the inner product of u and v modulo f
ZZX inner_prod(const vector<ZZX>& u, const vector<ZZX>& v, const ZZX& f){
	ZZX s = conv<ZZX>(0);
	ZZX result = conv<ZZX>(0);
	long l = u.size();
	if (u.size() != v.size()){
		cout << "ERROR: inner_prod: u.size() = " << u.size()
			 << " but v.size() = " << v.size() << endl;
		exit(1);
	}
	for (long i = 0; i < l; i++){
		// result += (u[i] * v[i] % f)
		MulMod(s, u[i], v[i], f);
		result += s;
	}
	return result;
}

// return the inner product of u and v modulo f
ZZX inner_prod(const vec_ZZ& u, const vector<ZZX>& v){
	ZZX s = conv<ZZX>(0);
	ZZX result = conv<ZZX>(0);
	long l = v.size();
	for (long i = 0; i < l; i++){
		result += u[i] * v[i];
	}
	return result;
}

// return the inner product of u and v modulo f
ZZ_pX inner_prod(const vector<ZZ_pX>& u, const vector<ZZ_pX>& v, const ZZ_pXModulus& f){
	ZZ_pX s = conv<ZZ_pX>(0);
	ZZ_pX result = conv<ZZ_pX>(0);
	long l = u.size();
	for (long i = 0; i < l; i++){
		// result += (u[i] * v[i] % f)
		MulMod(s, u[i], v[i], f);
		result += s;
	}
	return result;
}



bool rand_bit(){
	return (rand() % 2);
}

/**
 * Return  round(a/n), that is, interpret a/n as an rational number
 * then return the closest integer to it.
 */
ZZ rounded_division(const ZZ& a, const ZZ& n){
	long signal = (a >= 0 ? 1 : -1);
	ZZ _a = a*signal;
	// interpret a = q*n + r with 0 <= r < n
	ZZ q = _a / n;
	ZZ r = _a % n;
	if (2*r < n)
		return signal * q;
	else
		return signal * (q + 1);
}

/**
 * Return ceil(a/n), that is, interpret a/n as an rational number
 * then return the minimun integer bigger than a/n.
 */
ZZ ceil_division(const ZZ& a, const ZZ& n){
	long signal = (a >= 0 ? 1 : -1);
	ZZ _a = a*signal;
	// interpret a = q*n + r with 0 <= r < n
	ZZ q = _a / n;
	ZZ r = _a % n;
	if (0 != r){
		if (signal > 0)
			return (q + 1);
		else
			return -q;
	}
	return signal * q;
}
/** Let f = sum_{i=0}^N f_i * x^i be a polynomial of degree n.
 *  This function returns
 *      sum_{i=0}^N round(f_i / n) * x^i
 * that is, it applies rounded_division(ZZ, ZZ) to each coefficient of f.
 */
ZZX rounded_division(const ZZX& f, const ZZ& n){
	ZZX result;
	ZZ f_i;
	for (long i = deg(f); i >= 0; i--){
		f_i = coeff(f, i);
		SetCoeff(result, i, rounded_division(f_i, n));
	}
	return result;
}


ZZ symmetric_mod(const ZZ& a, const ZZ& n){ // returns a % n using the set ]-c/2, c/2] as Z/nZ
	ZZ b(a % n);
	if(2*b > n)
		b = b - n;
	return b;
}

/**
 *	 Reduces the real number x modulo n, that is, returns a real value r such
 * that -n/2 < r <= n/2 and a - r is an integer that is a multiple of n.
 */
RR mod_on_R(const RR& x, const ZZ& n){
	RR n_real(conv<RR>(n));
	ZZ int_part(TruncToZZ(x));
	RR decimal_part(x - conv<RR>(int_part));
	// reduces the integer part to ]-n/2,...,n/2]
	int_part = int_part % n;
	if(2*int_part > n)
		int_part = int_part - n;

	RR mod_n = conv<RR>(int_part) + decimal_part;

	if (2*mod_n < -n_real) // if reduced value is smaller than n/2
		mod_n += n_real;  // then add n
	
	if (2*mod_n > n_real)
		mod_n -= n_real;

	return mod_n;
}


/**    Set words_vec to the vector (x0, x1, ..., x_{l-1})
 * representing the decomposition of x in base b.
 */
void decompose_ZZ(vec_ZZ& words_vec, const ZZ& _x, long l, ZZ b){
	ZZ x(_x);
	long sign = (x < 0 ? -1 : 1);
	x *= sign;
	for (int j = 0; j < l; j++){
		words_vec[j] = sign * (x % b);
		x /= b;
	}
}

/**    Set words_vec to the vector (y0, y1, ..., y_{l-1})
 * representing the decomposition of f in base b, that is, each yi is a
 * polynomial with degree smaller than N and coefficients in ]-b, b[, and
 * the sum y0*b^0 + y1*b^1 + ... + y_{l-1}*b^{l-1} equals f
 */
void decompose_ZZX(vector<ZZX>& words_vec, const ZZX& f, long N, long l, ZZ b){
	ZZX pow_x(1);
	ZZX x; SetX(x);
	vec_ZZ decomp_coef; decomp_coef.SetLength(l);
	// decomp_coef = g^-1(f_0)
	decompose_ZZ(decomp_coef, f[0], l, b);
	for (long j = 0; j < l; j++)
		words_vec[j] = conv<ZZX>(decomp_coef[j]);

	for(long i = 1; i < N; i++){
		pow_x *= x; // pow_x = x^i 
		if (coeff(f, i) != 0){
			decompose_ZZ(decomp_coef, f[i], l, b);
			// After this loop: words_vec += x^i * g^-1(f_i)
			for (long j = 0; j < l; j++){
				words_vec[j] += pow_x * decomp_coef[j];
			}
		}
	}
}


/**    Set words_vec to the vector (y0, y1, ..., y_{l-1})
 * representing the decomposition of f in base b, that is, each yi is a
 * polynomial with degree smaller than N and coefficients in ]-b, b[, and
 * the sum y0*b^0 + y1*b^1 + ... + y_{l-1}*b^{l-1} equals f
 */
void decompose_ZZ_pX(vector<ZZ_pX>& words_vec, const ZZ_pX& f, long N, long l, ZZ b){
	ZZ_pX pow_x(1);
	ZZ_pX x; SetX(x);
	vec_ZZ decomp_coef; decomp_coef.SetLength(l);
	// decomp_coef = g^-1(f_0)
	decompose_ZZ(decomp_coef, conv<ZZ>(f[0]), l, b);
	for (long j = 0; j < l; j++)
		words_vec[j] = conv<ZZ_p>(decomp_coef[j]);

	for(long i = 1; i < N; i++){
		pow_x *= x; // pow_x = x^i 
		if (coeff(f, i) != 0){
			decompose_ZZ(decomp_coef, conv<ZZ>(f[i]), l, b);
			// After this loop: words_vec += x^i * g^-1(f_i)
			for (long j = 0; j < l; j++){
				words_vec[j] += pow_x * conv<ZZ_p>(decomp_coef[j]);
			}
		}
	}
}

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
vec_ZZ multiply_by_G(const vec_ZZ& v, long l, ZZ b){
	long N = v.length();
	vec_ZZ Gv; Gv.SetLength(N * l);
	for(long i = 0; i < N; i++){
		ZZ pow_b = ZZ(1);
		for(long j = 0; j < l; j++){
			Gv[i*l + j] = pow_b * v[i];
			pow_b *= b;
		}
	}
	return Gv;
}


/*
 *      Set words_vec to G^-1(vec), that is, assuming that vec is an n-dimensional vector
 * and words_vec an (n*l)-dimensional, then words_vec[0:l-1] corresponds to the base-b decomposition
 * vec[0], bits_vec[l:2*l-1] corresponds to the base-b decomposition of vec[1], and so on.
 * 		Thus l must have a value such that max(abs(vec[i])) < b^l.
 */
void invG(vec_ZZ& words_vec, const vec_ZZ& vec, long l, ZZ b){
	long n = vec.length();
	for(int i = 0; i < n; i++){
		ZZ x = vec[i];
		long sign = (x < 0 ? -1 : 1);
		x *= sign;
		for (int j = 0; j < l; j++){
			words_vec[i*l + j] = sign * (x % b);
			x /= b;
		}
	}
}

void invG_mat(mat_ZZ& words_mat, mat_ZZ& mat, long l, ZZ b){
	long nrows = mat.NumRows();
	for (long i = 0; i < nrows; i++){
		invG(words_mat[i], mat[i], l, b);
	}
}

vec_ZZ coefficient_vector(const ZZX& f, long N){
	vec_ZZ phi_f; phi_f.SetLength(N);
	for(long i = 0; i < N; i++){
		phi_f[i] = coeff(f, i);
	}
	return phi_f;
}

/* Receives a polynomial fmod of degree N and a polynomial f
 * of degree smaller than N.
 * Returns an NxN matrix F such that, for any polynomial g, the following holds:
 * 		 phi(g * f % fmod) = phi(g) * F
 * where phi is the function coefficient_vector
 **/
mat_ZZ circulant_matrix(const ZZX& f, const ZZX& fmod){
	long N = deg(fmod);
	mat_ZZ F; F.SetDims(N, N);
	vec_ZZ coefs; 
	ZZX h = f;
	for (long i = 0; i < N; i++){
		coefs = coefficient_vector(h, N);
		for (long j = 0; j < N; j++){
			F[i][j] = coeff(h, j);
		}
		h = MulByXMod(h, fmod); // h = h * x % fmod
	}
	return F;
}


/* --------------------------------------------------------------------
* Auxiliary functions used by several encryption schemes based on AGCD
* --------------------------------------------------------------------- */

ZZ sample_r(long rho){
	ZZ ri;
	RandomBits(ri, rho-1);
	if (rand_bit())
		return ri;
	else
		return -1 * ri;
}

ZZX sample_rx(long rho, long N){
	ZZX r;
	for (long i = N-1; i >= 0; i--){
		SetCoeff(r, i, sample_r(rho));
	}
	return r;
}

vector<ZZX> sample_vec_rx(long rho, long N, long l){
	vector<ZZX> v(l);
	for (long i = 0; i < l; i++)
		v[i] = sample_rx(rho, N);
	return v;
}

ZZ sample_q(long gamma, long eta){
	ZZ q = RandomBits_ZZ(gamma - eta);
	return q;
}

ZZX sample_qx(long gamma, long eta, long N){
	ZZX q;
	for (long i = N-1; i >= 0; i--){
		SetCoeff(q, i, sample_q(gamma, eta));
	}
	return q;
}

vector<ZZX> sample_vec_qx(long gamma, long eta, long N, long l){
	vector<ZZX> v(l);
	for (long i = 0; i < l; i++)
		v[i] = sample_qx(gamma, eta, N);
	return v;
}

vector<ZZX>& operator+= (vector<ZZX>& vec, const vector<ZZX>& u){
	for (long i = 0; i < vec.size(); i++)
		vec[i] += u[i];
	return vec;	
}
vector<ZZX>& operator*= (vector<ZZX>& vec, const ZZ& t){
	for (long i = 0; i < vec.size(); i++)
		vec[i] *= t;
	return vec;	
}


/**  Returns true if p is (probably) prime. 
* This function is essentially the same as the one in in NTL's documentation:
* http://phd-sid.ethz.ch/debian/ntl/ntl-6.0.0/doc/tour-ex1.html
**/
bool is_prime(const ZZ& n) {
   if (n <= 1) return 0;

   // first, perform trial division by primes up to 2000
   PrimeSeq s;  // a class for quickly generating primes in sequence
   long p;
   p = s.next();  // first prime is always 2
   while (p && p < 2000) {
      if ((n % p) == 0) return (n == p);
      p = s.next();  
   }

   // second, perform t Miller-Rabin tests
   long t = 12;
   ZZ x;
   for (long i = 0; i < t; i++) {
      x = RandomBnd(n); // random number between 0 and n-1
      if (MillerWitness(n, x))
         return false;
   }
   return true;
}


NTL::ZZ find_modulus(int max0, int max1, long int N){
	long int bit_len_k = max0 + max1;
	ZZ k = NextPrime(ZZ(1) << bit_len_k); 
	ZZ p = N*k + 1;
	while (! is_prime(p)){
		k = NextPrime(k+2);
		p = N*k + 1;
	}
	return p;
}

NTL::ZZ find_modulus(int max0, int max1, long int N, long int L){
	return find_modulus(max0 + L, max1, N);
}

bool is_primitive_root(NTL::ZZ g, long int N, NTL::ZZ p) {
	if (g <= 1)
		return false;

	if (PowerMod(g, ZZ(N), p) != 1)
		return false;

	// test if order(g) is some power of two smaller than N
	ZZ pow_g = g*g % p;
	for (long int Ni = 2; Ni < N; Ni *= 2) { // each iteration tests Ni = 2^i
		if (pow_g == 1) // test if order(g) = Ni
			return false;
		pow_g = pow_g * pow_g % p;
	}
	return true;
}

NTL::ZZ find_Nth_primitive_root(long int N, NTL::ZZ p){
	ZZ k = (p - 1) / N;
	ZZ g = ZZ(2);
	while(g < p){
		if (is_primitive_root(g, N, p))
			return g;

		ZZ g_k = PowerMod(g, k, p);
		if (is_primitive_root(g_k, N, p))
			return g_k;
		// otherwise, order(g) = 2^i * k for i < log(N), thus, test next g
		g++;
	}
	return g;
}

long reverse_bits(long int x, long int logN) {
	long int result = 0;
	for (long int i = 0; i < logN; x/= 2, i++)
		result = (result << 1) | (x & 1U);
	return result;
}

void coordinatewise_prod(const vec_ZZ u, const vec_ZZ v, vec_ZZ& output){
	for (long int i = 0; i < output.length(); i++)
		output[i] = (u[i] * v[i]);
}

void coordinatewise_prod(const vec_ZZ u, const vec_ZZ v, vec_ZZ& output, const ZZ& p){
	for (long int i = 0; i < output.length(); i++)
		output[i] = (u[i] * v[i]) % p;
}


void apply_powers(NTL::vec_ZZ& u, const NTL::ZZ& b, const NTL::ZZ& p){
	ZZ pow_b(1);
	for (long int i = 0; i < u.length(); i++){
		u[i] = u[i] * pow_b % p;
		pow_b = pow_b * b % p;
	}
}

void centralized_mod(NTL::vec_ZZ& u, const NTL::ZZ& p){
	for (long int i = 0; i < u.length(); i++){
		u[i] = u[i] % p;
		if (u[i] * 2 > p)
			u[i] -= p;
	}
}


NTL::ZZX vec_ZZ_to_poly(const NTL::vec_ZZ& u){
	ZZX f;
	for (long int i = 0; i < u.length(); i++) 
		SetCoeff(f, i, u[i]);
	return f;
}

/* ------------------------------------------------------
 * 			Random values generation
 * -----------------------------------------------------*/
vec_ZZ random_vec(int dimension, long int bound){
	vec_ZZ vec;
	vec.SetLength(dimension);
	for (int i = 0; i < dimension; i++){
		vec[i] = rand() % bound; // random number between 0 and bound-1
		if (rand() % 2)
			vec[i] *= -1;
	}
	return vec;
}

ZZX random_poly(long int degree, long int bitlen){
	ZZX a;
	for (long int i = 0; i <= degree; i++)
		SetCoeff(a, i, RandomBits_ZZ(bitlen));
	return a;
}

vector<ZZX> random_vec_of_polys(long int m, int d, long int bitlen){
	vector<ZZX> vec = vector<ZZX>(d);
	for (int i = 0; i < d; i++){
		vec[i] = random_poly(m, bitlen);
	}
	return vec;
}

