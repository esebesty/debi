
// The fisher exact test is taken from RMath, Mathlib : A C Library of Special Functions

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998      Ross Ihaka
 *  Copyright (C) 2004-2009 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
*/

namespace Rmath{

#ifndef M_2PI
#define M_2PI           6.283185307179586476925286766559        /* 2*pi */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI   0.918938533204672741780329736406        /* log(sqrt(2*pi)) */
#endif

inline long double lgammafn(double x){return ::lgamma(x);}
 double stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const double sferr_halves[31] = {
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,  /* 0.5 */
        0.0810614667953272582196702,  /* 1.0 */
        0.0548141210519176538961390,  /* 1.5 */
        0.0413406959554092940938221,  /* 2.0 */
        0.03316287351993628748511048, /* 2.5 */
        0.02767792568499833914878929, /* 3.0 */
        0.02374616365629749597132920, /* 3.5 */
        0.02079067210376509311152277, /* 4.0 */
        0.01848845053267318523077934, /* 4.5 */
        0.01664469118982119216319487, /* 5.0 */
        0.01513497322191737887351255, /* 5.5 */
        0.01387612882307074799874573, /* 6.0 */
        0.01281046524292022692424986, /* 6.5 */
        0.01189670994589177009505572, /* 7.0 */
        0.01110455975820691732662991, /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
        nn = n + n;
        if (nn == (int)nn) return(sferr_halves[(int)nn]);
        return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double bd0(double x, double np)
{
    double ej, s, s1, v;
    int j;

    if (fabs(x-np) < 0.1*(x+np)) {
        v = (x-np)/(x+np);
        s = (x-np)*v;/* s using v -- change by MM */
        ej = 2*x*v;
        v = v*v;
        for (j=1; ; j++) { /* Taylor series */
            ej *= v;
            s1 = s+ej/((j<<1)+1);
            if (s1==s) /* last term was effectively 0 */
                return(s1);
            s = s1;
        }
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}


double lfastchoose(double n, double k)
{
	return lgammafn(n + 1.0) - lgammafn(k + 1.0) - lgammafn(n - k + 1.0);
}

void ml_error(int n)
{
    switch(n) {

    case ME_NONE:
        errno = 0;
        break;

    case ME_DOMAIN:
      throw std::domain_error("Bmath domain error");
      break;
    case ME_NOCONV:
      throw std::domain_error("failed to converge");
      break;

    case ME_RANGE:
      throw std::range_error("Bmath range error");
      break;

    default:
      throw std::logic_error("call to Bmath::ml_error with unknown error");
      break;
    }
}

double dbinom_raw(double x, double n, double p, double q, int give_log)
{
    double lf, lc;

    if (p == 0) return((x == 0) ? R_D__1 : R_D__0);
    if (q == 0) return((x == n) ? R_D__1 : R_D__0);

    if (x == 0) {
	if(n == 0) return R_D__1;
	lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
	return( R_D_exp(lc) );
    }
    if (x == n) {
	lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
	return( R_D_exp(lc) );
    }
    if (x < 0 || x > n) return( R_D__0 );

    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    lf = log(M_2PI) + log(x) + log1p(- x/n);

    return R_D_exp(lc - 0.5*lf);
}


double dhyper(double x, double r, double b, double n, int give_log)
{
    double p, q, p1, p2, p3;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(r) || ISNAN(b) || ISNAN(n))
	return x + r + b + n;
#endif

    if (R_D_negInonint(r) || R_D_negInonint(b) || R_D_negInonint(n) || n > r+b)
	ML_ERR_return_NAN;
    if (R_D_negInonint(x))
	return(R_D__0);

    x = R_D_forceint(x);
    r = R_D_forceint(r);
    b = R_D_forceint(b);
    n = R_D_forceint(n);

    if (n < x || r < x || n - x > b) return(R_D__0);
    if (n == 0) return((x == 0) ? R_D__1 : R_D__0);

    p = ((double)n)/((double)(r+b));
    q = ((double)(r+b-n))/((double)(r+b));

    p1 = dbinom_raw(x,	r, p,q,give_log);
    p2 = dbinom_raw(n-x,b, p,q,give_log);
    p3 = dbinom_raw(n,r+b, p,q,give_log);

    return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
}

static double pdhyper (double x, double NR, double NB, double n, int log_p)
{
/*
 * Calculate
 *
 *	    phyper (x, NR, NB, n, TRUE, FALSE)
 *   [log]  ----------------------------------
 *	       dhyper (x, NR, NB, n, FALSE)
 *
 * without actually calling phyper.  This assumes that
 *
 *     x * (NR + NB) <= n * NR
 *
 */
    long double sum = 0;
    long double term = 1;

    while (x > 0 && term >= DBL_EPSILON * sum) {
	term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
	sum += term;
	x--;
    }

    return log_p ? log1p(sum) : 1 + sum;
}


/* FIXME: The old phyper() code was basically used in ./qhyper.c as well
 * -----  We need to sync this again!
*/
double phyper (double x, double NR, double NB, double n,
	       int lower_tail, int log_p)
{
/* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */

    double d, pd;

#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(NR) || ISNAN(NB) || ISNAN(n))
	return x + NR + NB + n;
#endif

    x = floor (x + 1e-7);
    NR = R_D_forceint(NR);
    NB = R_D_forceint(NB);
    n  = R_D_forceint(n);

    if (NR < 0 || NB < 0 || !R_FINITE(NR + NB) || n < 0 || n > NR + NB)
	ML_ERR_return_NAN;

    if (x * (NR + NB) > n * NR) {
	/* Swap tails.	*/
	double oldNB = NB;
	NB = NR;
	NR = oldNB;
	x = n - x - 1;
	lower_tail = !lower_tail;
    }

    if (x < 0)
	return R_DT_0;
    if (x >= NR || x >= n)
	return R_DT_1;

    d  = dhyper (x, NR, NB, n, log_p);
    pd = pdhyper(x, NR, NB, n, log_p);

    return log_p ? R_DT_Log(d + pd) : R_D_Lval(d * pd);
}
}
//////////////////////////////////////////////////////

