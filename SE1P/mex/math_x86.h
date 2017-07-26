#ifndef __MATH_X86_H_
#define __MATH_X86_H_

#include <immintrin.h>
#include <math.h>

#ifdef __AVX__
static inline __m256d
se_mm256_inv_pd(__m256d x)
{
    const __m256d two  = _mm256_set1_pd(2.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    __m256d lu = _mm256_cvtps_pd(_mm_rcp_ps( _mm256_cvtpd_ps(x)));

    /* Perform two N-R steps for double precision */
    lu         = _mm256_mul_pd(lu, _mm256_sub_pd(two, _mm256_mul_pd(x, lu)));
    return _mm256_mul_pd(lu, _mm256_sub_pd(two, _mm256_mul_pd(x, lu)));
}

/* 1.0/x, 256-bit wide */
static inline __m256
se_mm256_inv_ps(__m256 x)
{
    const __m256 two = _mm256_set1_ps(2.0f);

    __m256       lu = _mm256_rcp_ps(x);

    return _mm256_mul_ps(lu, _mm256_sub_ps(two, _mm256_mul_ps(lu, x)));
}

static inline __m256d
se_mm256_abs_pd(__m256d x)
{
  const __m256d signmask  = _mm256_castsi256_pd( _mm256_set_epi32(0x7FFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF,
								  0x7FFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF) );
  
  return _mm256_and_pd(x, signmask);
}

static inline __m256
se_mm256_abs_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );

    return _mm256_and_ps(x, signmask);
}


#ifdef SINGLE
static __m256
se_mm256_log_ps(__m256 x)
{
    const __m256  expmask    = _mm256_castsi256_ps( _mm256_set1_epi32(0x7F800000) );
    const __m128i expbase_m1 = _mm_set1_epi32(127-1); /* We want non-IEEE format */
    const __m256  half       = _mm256_set1_ps(0.5f);
    const __m256  one        = _mm256_set1_ps(1.0f);
    const __m256  invsq2     = _mm256_set1_ps(1.0f/sqrt(2.0f));
    const __m256  corr1      = _mm256_set1_ps(-2.12194440e-4f);
    const __m256  corr2      = _mm256_set1_ps(0.693359375f);

    const __m256  CA_1        = _mm256_set1_ps(0.070376836292f);
    const __m256  CB_0        = _mm256_set1_ps(1.6714950086782716f);
    const __m256  CB_1        = _mm256_set1_ps(-2.452088066061482f);
    const __m256  CC_0        = _mm256_set1_ps(1.5220770854701728f);
    const __m256  CC_1        = _mm256_set1_ps(-1.3422238433233642f);
    const __m256  CD_0        = _mm256_set1_ps(1.386218787509749f);
    const __m256  CD_1        = _mm256_set1_ps(0.35075468953796346f);
    const __m256  CE_0        = _mm256_set1_ps(1.3429983063133937f);
    const __m256  CE_1        = _mm256_set1_ps(1.807420826584643f);

    __m256        fexp;
    __m256i       iexp;
    __m128i       iexp128a, iexp128b;
    __m256        mask;
    __m256i       imask;
    __m128i       imask128a, imask128b;
    __m256        x2;
    __m256        y;
    __m256        pA, pB, pC, pD, pE, tB, tC, tD, tE;

    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp  = _mm256_and_ps(x, expmask);
    iexp  = _mm256_castps_si256(fexp);

    iexp128b = _mm256_extractf128_si256(iexp, 0x1);
    iexp128a = _mm256_castsi256_si128(iexp);

    iexp128a  = _mm_srli_epi32(iexp128a, 23);
    iexp128b  = _mm_srli_epi32(iexp128b, 23);
    iexp128a  = _mm_sub_epi32(iexp128a, expbase_m1);
    iexp128b  = _mm_sub_epi32(iexp128b, expbase_m1);

    x     = _mm256_andnot_ps(expmask, x);
    x     = _mm256_or_ps(x, one);
    x     = _mm256_mul_ps(x, half);

    mask  = _mm256_cmp_ps(x, invsq2, _CMP_LT_OQ);

    x     = _mm256_add_ps(x, _mm256_and_ps(mask, x));
    x     = _mm256_sub_ps(x, one);

    imask = _mm256_castps_si256(mask);

    imask128b = _mm256_extractf128_si256(imask, 0x1);
    imask128a = _mm256_castsi256_si128(imask);

    iexp128a  = _mm_add_epi32(iexp128a, imask128a);
    iexp128b  = _mm_add_epi32(iexp128b, imask128b);

    iexp  = _mm256_castsi128_si256(iexp128a);
    iexp  = _mm256_insertf128_si256(iexp, iexp128b, 0x1);

    x2    = _mm256_mul_ps(x, x);

    pA    = _mm256_mul_ps(CA_1, x);
    pB    = _mm256_mul_ps(CB_1, x);
    pC    = _mm256_mul_ps(CC_1, x);
    pD    = _mm256_mul_ps(CD_1, x);
    pE    = _mm256_mul_ps(CE_1, x);
    tB    = _mm256_add_ps(CB_0, x2);
    tC    = _mm256_add_ps(CC_0, x2);
    tD    = _mm256_add_ps(CD_0, x2);
    tE    = _mm256_add_ps(CE_0, x2);
    pB    = _mm256_add_ps(pB, tB);
    pC    = _mm256_add_ps(pC, tC);
    pD    = _mm256_add_ps(pD, tD);
    pE    = _mm256_add_ps(pE, tE);

    pA    = _mm256_mul_ps(pA, pB);
    pC    = _mm256_mul_ps(pC, pD);
    pE    = _mm256_mul_ps(pE, x2);
    pA    = _mm256_mul_ps(pA, pC);
    y     = _mm256_mul_ps(pA, pE);

    fexp  = _mm256_cvtepi32_ps(iexp);
    y     = _mm256_add_ps(y, _mm256_mul_ps(fexp, corr1));

    y     = _mm256_sub_ps(y, _mm256_mul_ps(half, x2));
    x2    = _mm256_add_ps(x, y);

    x2    = _mm256_add_ps(x2, _mm256_mul_ps(fexp, corr2));

    return x2;
}

#endif //SINGLE


#elif __SSE4_2__
static inline __m128d
se_mm_inv_pd(__m128d x)
{
    const __m128d two  = _mm_set1_pd(2.0);

    // Lookup instruction only exists in single precision, convert back and forth... //
    __m128d lu = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));

    // Perform two N-R steps for double precision //
    lu         = _mm_mul_pd(lu, _mm_sub_pd(two, _mm_mul_pd(x, lu)));
    return _mm_mul_pd(lu, _mm_sub_pd(two, _mm_mul_pd(x, lu)));
}

static inline __m128d
se_mm_abs_pd(__m128d x)
{
  const __m128d signmask  = _mm_castsi128_pd( _mm_set_epi32(0x7FFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF) );

    return _mm_and_pd(x, signmask);
}

static __m128d
se_mm_exp_pd(__m128d exparg)
{
    const __m128d argscale = _mm_set1_pd(1.4426950408889634073599);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128d arglimit = _mm_set1_pd(1022.0);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m128d invargscale0  = _mm_set1_pd(6.93145751953125e-1);
    const __m128d invargscale1  = _mm_set1_pd(1.42860682030941723212e-6);

    const __m128d P2       = _mm_set1_pd(1.26177193074810590878e-4);
    const __m128d P1       = _mm_set1_pd(3.02994407707441961300e-2);
    /* P0 == 1.0 */
    const __m128d Q3       = _mm_set1_pd(3.00198505138664455042E-6);
    const __m128d Q2       = _mm_set1_pd(2.52448340349684104192E-3);
    const __m128d Q1       = _mm_set1_pd(2.27265548208155028766E-1);
    /* Q0 == 2.0 */
    const __m128d one      = _mm_set1_pd(1.0);
    const __m128d two      = _mm_set1_pd(2.0);

    __m128d       valuemask;
    __m128i       iexppart;
    __m128d       fexppart;
    __m128d       intpart;
    __m128d       x, z, z2;
    __m128d       PolyP, PolyQ;

    x             = _mm_mul_pd(exparg, argscale);

    iexppart  = _mm_cvtpd_epi32(x);
    intpart   = _mm_round_pd(x, _MM_FROUND_TO_NEAREST_INT);

    /* The two lowest elements of iexppart now contains 32-bit numbers with a correctly biased exponent.
     * To be able to shift it into the exponent for a double precision number we first need to
     * shuffle so that the lower half contains the first element, and the upper half the second.
     * This should really be done as a zero-extension, but since the next instructions will shift
     * the registers left by 52 bits it doesn't matter what we put there - it will be shifted out.
     * (thus we just use element 2 from iexppart).
     */
    iexppart  = _mm_shuffle_epi32(iexppart, _MM_SHUFFLE(2, 1, 2, 0));

    /* Do the shift operation on the 64-bit registers */
    iexppart  = _mm_add_epi32(iexppart, expbase);
    iexppart  = _mm_slli_epi64(iexppart, 52);

    valuemask = _mm_cmpge_pd(arglimit, se_mm_abs_pd(x));
    fexppart  = _mm_and_pd(valuemask, _mm_castsi128_pd(iexppart));

    z         = _mm_sub_pd(exparg, _mm_mul_pd(invargscale0, intpart));
    z         = _mm_sub_pd(z, _mm_mul_pd(invargscale1, intpart));

    z2        = _mm_mul_pd(z, z);

    PolyQ     = _mm_mul_pd(Q3, z2);
    PolyQ     = _mm_add_pd(PolyQ, Q2);
    PolyP     = _mm_mul_pd(P2, z2);
    PolyQ     = _mm_mul_pd(PolyQ, z2);
    PolyP     = _mm_add_pd(PolyP, P1);
    PolyQ     = _mm_add_pd(PolyQ, Q1);
    PolyP     = _mm_mul_pd(PolyP, z2);
    PolyQ     = _mm_mul_pd(PolyQ, z2);
    PolyP     = _mm_add_pd(PolyP, one);
    PolyQ     = _mm_add_pd(PolyQ, two);

    PolyP     = _mm_mul_pd(PolyP, z);

    z         = _mm_mul_pd(PolyP, se_mm_inv_pd(_mm_sub_pd(PolyQ, PolyP)));
    z         = _mm_add_pd(one, _mm_mul_pd(two, z));

    z         = _mm_mul_pd(z, fexppart);

    return z;
}

static __m256
se_mm256_exp_ps(__m256 exparg)
{
    const __m256  argscale      = _mm256_set1_ps(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m256  arglimit      = _mm256_set1_ps(126.0f);
    const __m128i expbase       = _mm_set1_epi32(127);

    const __m256  invargscale0  = _mm256_set1_ps(0.693359375f);
    const __m256  invargscale1  = _mm256_set1_ps(-2.12194440e-4f);

    const __m256  CE5           = _mm256_set1_ps(1.9875691500e-4f);
    const __m256  CE4           = _mm256_set1_ps(1.3981999507e-3f);
    const __m256  CE3           = _mm256_set1_ps(8.3334519073e-3f);
    const __m256  CE2           = _mm256_set1_ps(4.1665795894e-2f);
    const __m256  CE1           = _mm256_set1_ps(1.6666665459e-1f);
    const __m256  CE0           = _mm256_set1_ps(5.0000001201e-1f);
    const __m256  one           = _mm256_set1_ps(1.0f);

    __m256        exparg2, exp2arg;
    __m256        pE0, pE1;
    __m256        valuemask;
    __m256i       iexppart;
    __m128i       iexppart128a, iexppart128b;
    __m256        fexppart;
    __m256        intpart;

    exp2arg = _mm256_mul_ps(exparg, argscale);

    iexppart  = _mm256_cvtps_epi32(exp2arg);
    intpart   = _mm256_round_ps(exp2arg, _MM_FROUND_TO_NEAREST_INT);

    iexppart128b = _mm256_extractf128_si256(iexppart, 0x1);
    iexppart128a = _mm256_castsi256_si128(iexppart);

    iexppart128a = _mm_slli_epi32(_mm_add_epi32(iexppart128a, expbase), 23);
    iexppart128b = _mm_slli_epi32(_mm_add_epi32(iexppart128b, expbase), 23);

    iexppart  = _mm256_castsi128_si256(iexppart128a);
    iexppart  = _mm256_insertf128_si256(iexppart, iexppart128b, 0x1);
    valuemask = _mm256_cmp_ps(arglimit, se_mm256_abs_ps(exp2arg), _CMP_GE_OQ);
    fexppart  = _mm256_and_ps(valuemask, _mm256_castsi256_ps(iexppart));

    /* Extended precision arithmetics */
    exparg    = _mm256_sub_ps(exparg, _mm256_mul_ps(invargscale0, intpart));
    exparg    = _mm256_sub_ps(exparg, _mm256_mul_ps(invargscale1, intpart));

    exparg2   = _mm256_mul_ps(exparg, exparg);

    pE1       = _mm256_mul_ps(CE5, exparg2);
    pE0       = _mm256_mul_ps(CE4, exparg2);
    pE1       = _mm256_add_ps(pE1, CE3);
    pE0       = _mm256_add_ps(pE0, CE2);
    pE1       = _mm256_mul_ps(pE1, exparg2);
    pE0       = _mm256_mul_ps(pE0, exparg2);
    pE1       = _mm256_add_ps(pE1, CE1);
    pE0       = _mm256_add_ps(pE0, CE0);
    pE1       = _mm256_mul_ps(pE1, exparg);
    pE0       = _mm256_add_ps(pE0, pE1);
    pE0       = _mm256_mul_ps(pE0, exparg2);
    exparg    = _mm256_add_ps(exparg, one);
    exparg    = _mm256_add_ps(exparg, pE0);

    exparg    = _mm256_mul_ps(exparg, fexppart);

    return exparg;
}


static __m128d
se_mm_log_pd(__m128d x)
{
    /* Same algorithm as cephes library */
    const __m128d expmask    = _mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );

    const __m128i expbase_m1 = _mm_set1_epi32(1023-1); /* We want non-IEEE format */

    const __m128d half       = _mm_set1_pd(0.5);
    const __m128d one        = _mm_set1_pd(1.0);
    const __m128d two        = _mm_set1_pd(2.0);
    const __m128d invsq2     = _mm_set1_pd(1.0/sqrt(2.0));

    const __m128d corr1      = _mm_set1_pd(-2.121944400546905827679e-4);
    const __m128d corr2      = _mm_set1_pd(0.693359375);

    const __m128d P5         = _mm_set1_pd(1.01875663804580931796e-4);
    const __m128d P4         = _mm_set1_pd(4.97494994976747001425e-1);
    const __m128d P3         = _mm_set1_pd(4.70579119878881725854e0);
    const __m128d P2         = _mm_set1_pd(1.44989225341610930846e1);
    const __m128d P1         = _mm_set1_pd(1.79368678507819816313e1);
    const __m128d P0         = _mm_set1_pd(7.70838733755885391666e0);

    const __m128d Q4         = _mm_set1_pd(1.12873587189167450590e1);
    const __m128d Q3         = _mm_set1_pd(4.52279145837532221105e1);
    const __m128d Q2         = _mm_set1_pd(8.29875266912776603211e1);
    const __m128d Q1         = _mm_set1_pd(7.11544750618563894466e1);
    const __m128d Q0         = _mm_set1_pd(2.31251620126765340583e1);

    const __m128d R2         = _mm_set1_pd(-7.89580278884799154124e-1);
    const __m128d R1         = _mm_set1_pd(1.63866645699558079767e1);
    const __m128d R0         = _mm_set1_pd(-6.41409952958715622951e1);

    const __m128d S2         = _mm_set1_pd(-3.56722798256324312549E1);
    const __m128d S1         = _mm_set1_pd(3.12093766372244180303E2);
    const __m128d S0         = _mm_set1_pd(-7.69691943550460008604E2);

   __m128d       fexp;
    __m128i       iexp;

    __m128d       mask1, mask2;
    __m128d       corr, t1, t2, q;
    __m128d       zA, yA, xA, zB, yB, xB, z;
    __m128d       polyR, polyS;
    __m128d       polyP1, polyP2, polyQ1, polyQ2;

    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp   = _mm_and_pd(x, expmask);
    iexp   = _mm_castpd_si128(fexp);
    iexp   = _mm_srli_epi64(iexp, 52);
    iexp   = _mm_sub_epi32(iexp, expbase_m1);
    iexp   = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(1, 1, 2, 0) );
    fexp   = _mm_cvtepi32_pd(iexp);

    x      = _mm_andnot_pd(expmask, x);
    x      = _mm_or_pd(x, one);
    x      = _mm_mul_pd(x, half);

    mask1     = _mm_cmpgt_pd(se_mm_abs_pd(fexp), two);
    mask2     = _mm_cmplt_pd(x, invsq2);

    fexp   = _mm_sub_pd(fexp, _mm_and_pd(mask2, one));

    /* If mask1 is set ('A') */
    zA     = _mm_sub_pd(x, half);
    t1     = _mm_blendv_pd( zA, x, mask2 );
    zA     = _mm_sub_pd(t1, half);
    t2     = _mm_blendv_pd( x, zA, mask2 );
    yA     = _mm_mul_pd(half, _mm_add_pd(t2, one));

    xA     = _mm_mul_pd(zA, se_mm_inv_pd(yA));
    zA     = _mm_mul_pd(xA, xA);

    /* EVALUATE POLY */
    polyR  = _mm_mul_pd(R2, zA);
    polyR  = _mm_add_pd(polyR, R1);
    polyR  = _mm_mul_pd(polyR, zA);
    polyR  = _mm_add_pd(polyR, R0);

    polyS  = _mm_add_pd(zA, S2);
    polyS  = _mm_mul_pd(polyS, zA);
    polyS  = _mm_add_pd(polyS, S1);
    polyS  = _mm_mul_pd(polyS, zA);
    polyS  = _mm_add_pd(polyS, S0);

    q      = _mm_mul_pd(polyR, se_mm_inv_pd(polyS));
    zA     = _mm_mul_pd(_mm_mul_pd(xA, zA), q);

    zA     = _mm_add_pd(zA, _mm_mul_pd(corr1, fexp));
    zA     = _mm_add_pd(zA, xA);
    zA     = _mm_add_pd(zA, _mm_mul_pd(corr2, fexp));

    /* If mask1 is not set ('B') */
    corr   = _mm_and_pd(mask2, x);
    xB     = _mm_add_pd(x, corr);
    xB     = _mm_sub_pd(xB, one);
    zB     = _mm_mul_pd(xB, xB);

    polyP1 = _mm_mul_pd(P5, zB);
    polyP2 = _mm_mul_pd(P4, zB);
    polyP1 = _mm_add_pd(polyP1, P3);
    polyP2 = _mm_add_pd(polyP2, P2);
    polyP1 = _mm_mul_pd(polyP1, zB);
    polyP2 = _mm_mul_pd(polyP2, zB);
    polyP1 = _mm_add_pd(polyP1, P1);
    polyP2 = _mm_add_pd(polyP2, P0);
    polyP1 = _mm_mul_pd(polyP1, xB);
    polyP1 = _mm_add_pd(polyP1, polyP2);

    polyQ2 = _mm_mul_pd(Q4, zB);
    polyQ1 = _mm_add_pd(zB, Q3);
    polyQ2 = _mm_add_pd(polyQ2, Q2);
    polyQ1 = _mm_mul_pd(polyQ1, zB);
    polyQ2 = _mm_mul_pd(polyQ2, zB);
    polyQ1 = _mm_add_pd(polyQ1, Q1);
    polyQ2 = _mm_add_pd(polyQ2, Q0);
    polyQ1 = _mm_mul_pd(polyQ1, xB);
    polyQ1 = _mm_add_pd(polyQ1, polyQ2);

    fexp   = _mm_and_pd(fexp, _mm_cmpneq_pd(fexp, _mm_setzero_pd()));

    q      = _mm_mul_pd(polyP1, se_mm_inv_pd(polyQ1));
    yB     = _mm_mul_pd(_mm_mul_pd(xB, zB), q);

    yB     = _mm_add_pd(yB, _mm_mul_pd(corr1, fexp));
    yB     = _mm_sub_pd(yB, _mm_mul_pd(half, zB));
    zB     = _mm_add_pd(xB, yB);
    zB     = _mm_add_pd(zB, _mm_mul_pd(corr2, fexp));

    z      = _mm_blendv_pd( zB, zA, mask1 );

    return z;
}



#endif // AVX

#endif // _MATH_X86_H_
