/** @file algebra.h
 *  @brief 
 * 
 * PetscVector specific file for external/petscvector/algebra/algebra.h
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_PETSCVECTOR_ALGEBRA_H
#define	PASC_PETSCVECTOR_ALGEBRA_H

#include "general/algebra/algebra.h"

/* vector */
#include "external/petscvector/algebra/vector/generalvector.h"

/* matrix */
//#include "external/petscvector/algebra/matrix/generalmatrixrhs.h"
#include "external/petscvector/algebra/matrix/generalmatrix.h"
#include "external/petscvector/algebra/matrix/blockgraphsparse.h"
#include "external/petscvector/algebra/matrix/blockgraphsparseTR.h"

/* feasibleset */
//#include "external/petscvector/algebra/feasibleset/generalfeasibleset.h"
#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"
#include "external/petscvector/algebra/feasibleset/simplex_local.h"

/* graph */
#include "external/petscvector/algebra/graph/bgmgraph.h"
#include "external/petscvector/algebra/graph/bgmgraphgrid1D.h"
#include "external/petscvector/algebra/graph/bgmgraphgrid2D.h"
#include "external/petscvector/algebra/graph/bgmgraphgrid3D.h"

/* integration */
#include "external/petscvector/algebra/integration/entropyintegration.h"
#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"
//#include "external/petscvector/algebra/integration/entropyintegrationcuba.h"

/* fem */
#include "external/petscvector/algebra/fem/fem1Dsum.h"
#include "external/petscvector/algebra/fem/fem1Dhat.h"
#include "external/petscvector/algebra/fem/fem2Dsum.h"
#include "external/petscvector/algebra/fem/fem2Dhat.h"


#endif
