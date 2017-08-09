/** @file algebra.h
 *  @brief 
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_ALGEBRA_H
#define	PASC_ALGEBRA_H

/* feasibleset */
#include "general/algebra/feasibleset/generalfeasibleset.h"
#include "general/algebra/feasibleset/simplex_lineqbound.h"
#include "general/algebra/feasibleset/simplex_local.h"

/* graph */
#include "general/algebra/graph/bgmgraph.h"
#include "general/algebra/graph/bgmgraphgrid1D.h"
#include "general/algebra/graph/bgmgraphgrid2D.h"
#include "general/algebra/graph/bgmgraphgrid3D.h"

/* integration */
#include "general/algebra/integration/entropyintegration.h"
#include "general/algebra/integration/entropyintegrationdlib.h"
#include "general/algebra/integration/entropyintegrationcuba.h"
#include "general/algebra/integration/entropyintegrationcudavegas.h"

/* matrix */
#include "general/algebra/matrix/generalmatrixrhs.h"
#include "general/algebra/matrix/generalmatrix.h"
#include "general/algebra/matrix/blockgraphsparse.h"
#include "general/algebra/matrix/blockgraphsparseTR.h"

/* vector */
#include "general/algebra/vector/generalvector.h"

/* fem */
#include "general/algebra/fem/fem.h"
#include "general/algebra/fem/fem1Dsum.h"
#include "general/algebra/fem/fem1Dhat.h"
#include "general/algebra/fem/fem2Dsum.h"
#include "general/algebra/fem/fem2Dhat.h"


#endif
