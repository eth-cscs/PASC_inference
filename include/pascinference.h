/*******************************************************************************
PASC INFERENCE library
Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
USI Lugano, 2016
lukas.pospisil@usi.ch

*******************************************************************************/

#ifndef PASCINFERENCE_H
#define	PASCINFERENCE_H


/* include common c++ header files */
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stack>
#include <limits>



/* include pascinference stuff */
#include "common.h"
#include "algebra.h"

#include "generaldata.h"
#include "generalresult.h"
#include "generalsolver.h"
#include "generalfeasibleset.h"



/*
#include "projection.h"
#include "qpsolver.h"
#include "theta.h"

#include "model.h"
#include "operations.h"
#include "inputoutput.h"
#include "data.h"
#include "gamma.h"
*/

/* Doxygen main page follows */

/*! \mainpage PASC HPC-Causality project
 *
 * Analysis of large amounts of economical data and data-driven inference of causality relations between different components of economical systems is one of the central problems in modern computational finance and economics. The task of proper mathematical description and adequate causality understanding for the economical data is hampered by the multiscale nature of the underlying processes, resulting from the presence of different temporal and spatial, i.e. regional, sectorial and global scales. Important questions thereby are: (i) an investigation of the mutual causality influences of different economic observables and their spatial (e.g., regional) and temporal (e.g., associated with the business cycle) evolution, (ii) identification of the most important exogenous impact factors that play a role in their dynamics, (iii) proper mathematical and statistical description of the influences coming from the unresolved/latent scales and factors.
 *
 * The main practical question addressed by the project will be to analyse one of the most comprehensive economical data bases available for credit risk analysis, containing the rating histories for 15'123 individual companies from 27-May-1977 to 8-Jan-2014 that were provided by Standard&Poor's. Implementation of data driven statistical methods for causality inference in such big data sets leads to a problem, which is intractable on single socket computing platforms. From the perspective of economic analysis, in this project we would like to identify whether there are some statistically significant causality relations between the credit rating migrations and, if yes, whether they do change in time or not (to a statistically and economically significant extent), e.g., due to the impact of either observable, or latent, economical factors, such as business cycle fluctuations, or market and funding liquidity conditions. Disentangling the effects of the causality network (contagion) from common risk factors (systematic risks) is of crucial importance for the design of regulations of financial markets, since these two forms of risk dependencies demand different types of supervisory activities. The project also aims at exploring how the inferred causality links are distributed in relation to characteristics of the companies, such as their industrial sector, financial size, etc. To answer the above economic questions, the project deploys data-driven causality inference methods. As another research goal, the results obtained with data-driven approaches will be compared with those of structural econometric methods based on large panel models with common stochastic factors and contagion effects. The trade-off between model flexibility and computational complexity will be investigated.
 *
 * Achieving the aims formulated in this project would extend the current prototypical and sequential methodology for data-driven causality inference. Developed algorithms and software would lay the groundwork for the future HPC-driven analysis of constantly growing sets of economical data and inference of very large causality networks subject to impacts of exogenous and latent economical factors beyond the limitations of currently available methods (e.g., beyond the current limitations of the Granger causality inference). Application of this new framework will provide an alternative network-driven multiscale perspective and new tools for understanding the interactions in very complex economical systems.
 *
 * 
 */




#endif


