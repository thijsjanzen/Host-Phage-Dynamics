/**************************   mersenne.cpp   **********************************
* Author:        Agner Fog
* Date created:  2001
* Last modified: 2008-11-16
* Project:       randomc.h
* Platform:      Any C++
* Description:
* Random Number generator of type 'Mersenne Twister'
*
* This random number generator is described in the article by
* M. Matsumoto & T. Nishimura, in:
* ACM Transactions on Modeling and Computer Simulation,
* vol. 8, no. 1, 1998, pp. 3-30.
* Details on the initialization scheme can be found at
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

#include "randomc.h"
#include <cmath>
#include <iostream>

CRandomMersenne rndgen(5); //the one random number generator

void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   int i, j, k;

   // Initialize
   Init0(19650218);

   if (NumSeeds <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= NumSeeds) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((double)(uint32_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_SUPPORTED
   // 64 bit integers available. Use multiply and shift method
   uint32_t interval;                    // Length of interval
   uint64_t longran;                     // Random bits * interval
   uint32_t iran;                        // Longran / 2^32
   uint32_t remainder;                   // Longran % 2^32

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRandom() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint32_t interval;                    // Length of interval
   uint32_t bran;                        // Random bits
   uint32_t iran;                        // bran / interval
   uint32_t remainder;                   // bran % interval

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32_t)0xFFFFFFFF / interval;
      if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int32_t)remainder + min;

#endif
}



void set_seed(int seed)
{
	rndgen.RandomInit(seed);
}

double uniform()
{
	return rndgen.Random();
}

int random_number(int n)
{
	return rndgen.IRandom(0,n-1);
}

double normal(double m, double s)
{
	return rndgen.normal(m,s);
}


double Expon(double lambda)
{
	return log(1 - rndgen.Random()) / (-1.0 * lambda);
}

double Binom(int n, double p)
{
	return rndgen.Binomial(n,p);
}





double CRandomMersenne::normal(double m, double s)
{
	double normal_x1;                   // first random coordinate (normal_x2 is member of class)
	double w;                           // radius
	
	if (normal_x2_valid) {              // we have a valid result from last call
		normal_x2_valid = 0;
		return normal_x2 * s + m;
	}
	
	// make two normally distributed variates by Box-Muller transformation
	do {
		normal_x1 = 2. * Random() - 1.;
		normal_x2 = 2. * Random() - 1.;
		w = normal_x1*normal_x1 + normal_x2*normal_x2;
	} while (w >= 1. || w < 1E-30);
	
	w = std::sqrt(std::log(w)*(-2./w));
	normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
	normal_x2_valid = 1;                // save normal_x2 for next call
	return normal_x1 * s + m;
}

/***********************************************************************
 Binomial distribution
 ***********************************************************************/
int CRandomMersenne::Binomial (int32_t n, double p) {
	/*
	 This function generates a random variate with the binomial distribution.
	 
	 Uses inversion by chop-down method for n*p < 35, and ratio-of-uniforms
	 method for n*p >= 35.
	 
	 For n*p < 1.E-6 numerical inaccuracy is avoided by poisson approximation.
	 */
	int inv = 0;                        // invert
	int32_t x;                          // result
	double np = n * p;
	
	if (p > 0.5) {                      // faster calculation by inversion
		p = 1. - p;  inv = 1;
	}
	if (n <= 0 || p <= 0) {
		if (n == 0 || p == 0) {
			return inv * n;  // only one possible result
		}
		// error exit
		std::cout << "Parameter out of range in binomial function";
		exit(1);
	}
	
	//------------------------------------------------------------------
	//                 choose method
	//------------------------------------------------------------------
	if (np < 35.) {
		if (np < 1.E-6) {
			// Poisson approximation for extremely low np
			x = PoissonLow(np);
		}
		else {
			// inversion method, using chop-down search from 0
			x = BinomialInver(n, p);
		}
	}
	else {
		// ratio of uniforms method
		x = BinomialRatioOfUniforms(n, p);
	}
	if (inv) {
		x = n - x;      // undo inversion
	}
	return x;
}

int32_t CRandomMersenne::PoissonLow(double L) {
	/*
	 This subfunction generates a random variate with the poisson
	 distribution for extremely low values of L.
	 
	 The method is a simple calculation of the probabilities of x = 1
	 and x = 2. Higher values are ignored.
	 
	 The reason for using this method is to avoid the numerical inaccuracies
	 in other methods.
	 */
	double d, r;
	d = sqrt(L);
	if (Random() >= d) return 0;
	r = Random() * d;
	if (r > L * (1.-L)) return 0;
	if (r > 0.5 * L*L * (1.-L)) return 1;
	return 2;
}


/***********************************************************************
 Subfunctions used by binomial
 ***********************************************************************/

int32_t CRandomMersenne::BinomialInver (int32_t n, double p) {
	/*
	 Subfunction for Binomial distribution. Assumes p < 0.5.
	 
	 Uses inversion method by search starting at 0.
	 
	 Gives overflow for n*p > 60.
	 
	 This method is fast when n*p is low.
	 */
	double f0, f, q;
	int32_t bound;
	double pn, r, rc;
	int32_t x, n1, i;
	
	// f(0) = probability of x=0 is (1-p)^n
	// fast calculation of (1-p)^n
	f0 = 1.;  pn = 1.-p;  n1 = n;
	while (n1) {
		if (n1 & 1) f0 *= pn;
		pn *= pn;  n1 >>= 1;
	}
	// calculate safety bound
	rc = (n + 1) * p;
	bound = (int32_t)(rc + 11.0*(sqrt(rc) + 1.0));
	if (bound > n) bound = n;
	q = p / (1. - p);
	
	while (1) {
		r = Random();
		// recursive calculation: f(x) = f(x-1) * (n-x+1)/x*p/(1-p)
		f = f0;  x = 0;  i = n;
		do {
			r -= f;
			if (r <= 0) return x;
			x++;
			f *= q * i;
			r *= x;       // it is faster to multiply r by x than dividing f by x
			i--;
		}
		while (x <= bound);
	}
}

double CRandomMersenne::LnFac(int32_t n) {
	// log factorial function. gives natural logarithm of n!
	
	// define constants
	static const double                 // coefficients in Stirling approximation
	C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
	C1 =  1./12.,
	C3 = -1./360.;
	static const int FAK_LEN = 1024;       // length of factorial table
	// C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
	// C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
	// static variables
	static double fac_table[FAK_LEN];   // table of ln(n!):
	static int initialized = 0;         // remember if fac_table has been initialized
	
	if (n < FAK_LEN) {
		if (n <= 1) {
			if (n < 0) {
				std::cout << "Parameter negative in LnFac function\n";
				exit(1);
			}
			return 0;
		}
		if (!initialized) {              // first time. Must initialize table
			// make table of ln(n!)
			double sum = fac_table[0] = 0.;
			for (int i=1; i<FAK_LEN; i++) {
				sum += log(double(i));
				fac_table[i] = sum;
			}
			initialized = 1;
		}
		return fac_table[n];
	}
	// not found in table. use Stirling approximation
	double  n1, r;
	n1 = n;  r  = 1. / n1;
	return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}


int32_t CRandomMersenne::BinomialRatioOfUniforms (int32_t n, double p) {
	/*
	 Subfunction for Binomial distribution. Assumes p < 0.5.
	 
	 Uses the Ratio-of-Uniforms rejection method.
	 
	 The computation time hardly depends on the parameters, except that it matters
	 a lot whether parameters are within the range where the LnFac function is
	 tabulated.
	 
	 Reference: E. Stadlober: "The ratio of uniforms approach for generating
	 discrete random variates". Journal of Computational and Applied Mathematics,
	 vol. 31, no. 1, 1990, pp. 181-189.
	 */
	double u;                           // uniform random
	double q1;                          // 1-p
	double np;                          // n*p
	double var;                         // variance
	double lf;                          // ln(f(x))
	double x;                           // real sample
	int32_t k;                          // integer sample
	
	static const double SHAT1 = 2.943035529371538573;    // 8/e
	static const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)
	
	if(bino_n_last != n || bino_p_last != p) {    // Set_up
		bino_n_last = n;
		bino_p_last = p;
		q1 = 1.0 - p;
		np = n * p;
		bino_mode = (int32_t)(np + p);             // mode
		bino_a = np + 0.5;                         // hat center
		bino_r1 = log(p / q1);
		bino_g = LnFac(bino_mode) + LnFac(n-bino_mode);
		var = np * q1;                             // variance
		bino_h = sqrt(SHAT1 * (var+0.5)) + SHAT2;  // hat width
		bino_bound = (int32_t)(bino_a + 6.0 * bino_h);// safety-bound
		if (bino_bound > n) bino_bound = n;        // safety-bound
	}
	
	while (1) {                                   // rejection loop
		u = Random();
		if (u == 0) continue;                      // avoid division by 0
		x = bino_a + bino_h * (Random() - 0.5) / u;
		if (x < 0. || x > bino_bound) continue;    // reject, avoid overflow
		k = (int32_t)x;                            // truncate
		lf = (k-bino_mode)*bino_r1+bino_g-LnFac(k)-LnFac(n-k);// ln(f(k))
		if (u * (4.0 - u) - 3.0 <= lf) break;      // lower squeeze accept
		if (u * (u - lf) > 1.0) continue;          // upper squeeze reject
		if (2.0 * log(u) <= lf) break;             // final acceptance
	}
	return k;
}








