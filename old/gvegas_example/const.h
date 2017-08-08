/* This is a -*- C -*- header. */
/* =================================================================== *\
   last updated : [2008/07/03] -11:59:39(JST)-
   by $Author: okamura $ on $Date: 2008/07/03 05:42:37 $ 
   $Name:  $-$Revision: 1.2 $-$State: Exp $
   $Source: /home/okamura/cvsroot/gMkMatrixE.uu.aa/const.h,v $ 
\* =================================================================== */


/*!
 * this file is for the mathematical and physical constant.
 */

#pragma once

/*-------- physical constants ---------*/
/*! Fermi CC */
//#define GF          (1.16637*1.0e-5)    /* GeV^-2 */
/*! hc */
#define hc          (0.197326968)       /* GeV.fm */
/*! Avogadro Num */
#define NA          (6.0221415*1.0e23)  /* 1/mol */
/*! cross section */
#define GeV2pb      (0.389379304*1.0e9)     /* GeV^-2 pbarn */
#define GeV2fb      (0.389379304*1.0e12)     /* GeV^-2 pbarn */

// alpha = 1/127
//#define hALPHA 7.87401575e-3
#define hALPHA 7.546772e-3

// ee2 = alpha * 4pi
//#define hEE2 9.89478e-2
#define hEE2 9.48356e-2
// ee = sqr(ee2)
//#define hEE 0.31460f
#define hEE 0.30795f

#define Clr 8

//---------------------------------------
//SU(3) coupling constant g^2/4pi = 0.118
//---------------------------------------
//#define GccCL   -1.2177158f   // left-handed gcc
//#define GccCR   -1.2177158f   // right-handed gcc
//#define GG       1.2177158f   // gcc at vvv

//---------------------------------------
//U(1) coupling constant g^2/4pi = 1/127
//---------------------------------------
//#define GccEL   -0.307953790484f   // left-handed gcc
//#define GccER   -0.307953790484f   // right-handed gcc


/*-------- mathematical constants ---------*/
/*! rad -> deg */
#define  Rad2Deg    57.2957795130823f
/*! deg -> rad */
#define  Deg2Rad    0.01745329251994f

/*! pi */
#define PI          3.14159265358979f
/*! pi/2 */
#define PIhalf      1.570796326794895f
/*! 2pi */
#define PI2         6.28318530717959f
#define TWOPI       6.28318530717959f
/*! pi^2 */
#define	PIPI        9.86960440108936f
/*! 1/pi */
#define	rPI         0.31830988618379f
/*! 1/2pi */
#define	rPI2        0.159154943091895f
#define	rTWOPI      0.159154943091895f
/*! sqrt(pi) */
#define	SqPI        1.77245385090552f
/*! 1/sqrt(pi) */
#define	rSqPI       0.56418958354776f

/*! sqrt(2) */
#define	SQRT2		1.41421356237310f
/*! 1/sqrt(2) */
#define	rSQRT2  	0.70710678118655f
/*! sqrt(3) */
#define	SQRT3		1.73205080756888f
/*! 1/sqrt(3) */
#define	rSQRT3  	0.57735026918963f
/*! sqrt(5) */
#define	SQRT5		2.23606797749979f
/*! 1/sqrt(5) */
#define	rSQRT5  	0.44721359549996f


/*! e */
#define	EXP         2.71828182845905f
/*! 1/e */
#define	rEXP        0.36787944117144f

/*! ln 10 */
#define	LN10		2.30258509299405f
/*! log  e */
#define	LOGE		0.43429448190325f

#define LOGHALF       -0.693147180559945f

#define zero           0.f

#define HALF           0.500000000000000f
#define THIRDS         0.333333333333333f
#define NINETH         0.111111111111111f
#define TWENTYSEVENTH  0.037037037037037f

