#ifndef __GAHE_PARAMS__
#define __GAHE_PARAMS___

#include <map>

class PreComputedParams {

	public:
	
	long lambda;
	long N;
	long eta;
	long rho;
	long gamma;
	long log_b;
	
	PreComputedParams(long lambda, long N, 
		long eta, long rho, long gamma, long log_b){
		this->lambda = lambda;
		this->N = N;
		this->eta = eta;
		this->rho = rho;
		this->gamma = gamma;
		this->log_b = log_b;
	}
};

std::map<long , PreComputedParams> params_100_private_x0 {
//{log_N, PreComputedParams(lambda,  N,   eta, rho, gamma, log_b)},
{  2, PreComputedParams(    100,     4,   100,  73,  2744,   7)},
{  3, PreComputedParams(    100,     8,   100,  73,  1372,   7)},
{ 4, PreComputedParams(     100,    16,   100,  73,   685,   7)},
{ 5, PreComputedParams(     100,    32,   100,  73,   343,   7)},
{ 6, PreComputedParams(     100,    64,   100,  72,   200,   11)},
{ 7, PreComputedParams(     100,   128,   100,  59,   200,   19)},
{ 8, PreComputedParams(     100,   256,   100,  42,   200,   36)},
{ 9, PreComputedParams(     100,   512,   100,  18,   200,   60)},
{10, PreComputedParams(     100,  1024,   100,   2,   200,   76)},
};


#endif
