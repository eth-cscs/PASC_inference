/*******************************************************************************
Conjugate Gradient Method for PASC INFERENCE library
Lukas Pospisil, 2016
lukas.pospisil@vsb.cz

* using MINLIN library (CSCS Lugano - Timothy Moroney, Ben Cumming)
* created during HPCCausality project (USI Lugano - Illia Horenko, Patrick Gagliardini, Will Sawyer)
*******************************************************************************/

#ifndef SOLVERCG_H
#define	SOLVERCG_H

#include <iostream>
#include "algebra/algebra.h"

namespace pascinference {

	/* unconstrained with initial approximation */
	template<class VectorType,class VectorTypeMatrix>
	VectorType cg(const GeneralMatrix<VectorTypeMatrix> &A, const VectorType &b, const VectorType &x0){

		VectorType x;
		x = x0; /* set approximation as initial */

		/* CG method */
		VectorType g(x0.size()); /* gradient */
		VectorType p(x0.size()); /* A-conjugate vector */
		VectorType Ap(x0.size()); /* A*p */

		int it = 0; /* iteration counter */
		int hess_mult = 0; /* number of hessian multiplications */
		double normg, alpha, beta, pAp, gg, gg_old;
	
		g = A*x; hess_mult += 1; g -= b; /* compute gradient */
		p = g; /* initial conjugate gradient */

		gg = dot(g,g);
		normg = std::sqrt(gg);

		while(normg > 0.001 && it < 10000){
			/* compute new approximation */

			Ap = A*p; hess_mult += 1;
			
			pAp = dot(Ap,p);
			alpha = gg/pAp; /* compute step-size */
			x -= alpha*p; /* set new approximation */

			/* compute gradient recursively */
			g -= alpha*Ap; 
			gg_old = gg;
			gg = dot(g,g);
			normg = std::sqrt(gg);
			
			/* compute new A-orthogonal vector */
			beta = gg/gg_old;
			p *= beta;
			p += g;
		
			std::cout << "it " << it << ": ||g|| = " << normg << std::endl;

			if(DEBUG_MODE >= 10){
				std::cout << "x = " << x << std::endl;
				std::cout << "g = " << g << std::endl;
				std::cout << "p = " << p << std::endl;
				std::cout << "Ap = " << Ap << std::endl;
				std::cout << "pAp = " << pAp << std::endl;
				std::cout << "alpha = " << alpha << std::endl;
				std::cout << "beta = " << beta << std::endl;
				std::cout << "gg = " << gg << std::endl;
				std::cout << "gg_old = " << gg_old << std::endl;
				std::cout << "normg = " << normg << std::endl;

				std::cout << "------------------------------------" << std::endl;
			}

				
			it += 1;

		}
		
		/* print output */
		std::cout << "------------------------" << std::endl;
		std::cout << " it_cg = " << it << std::endl;
		std::cout << " norm_g = " << normg << std::endl;
		std::cout << " hess_mult = " << hess_mult << std::endl;

		return x;
	}


} /* end namespace */




#endif
