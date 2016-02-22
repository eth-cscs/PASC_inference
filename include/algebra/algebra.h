/* these files includes all fun with vectors and matrices:
 * - define GlobalVector, HostVector, DeviceVector
 * - define GeneralMatrix and operations with it
 * - define projection onto Bx=1, x>=0	
*/

#include "common.h"
#include "algebra/generalmatrix.h"

// TODO: here I should deal with USE_MINLIN, USE_PETSC, USE_GPU etc.

typedef petscvector::PetscVector PetscVector;
typedef minlin::threx::HostVector<double> MinlinHostVector;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;


namespace pascinference {

	/* GlobalVector from PetscVector */
	class GlobalVector : public PetscVector {
		
		public:
			/* constructors */
			GlobalVector(int N): PetscVector(N) {}
			
			/* matrix-vector multiplication with General matrix */
			GlobalVector &operator=(GeneralMatrixRHS<GlobalVector> rhs){
				rhs.matmult(*this);
				return *this;
			}
			
	};




	/* HostVector from minlin::HostVector */
	class HostVector : public MinlinHostVector {
		
		public:
			HostVector(int N): MinlinHostVector(N) {} /* constructor */
			
	};



	/* DeviceVector from minlin::DeviceVector */
	class DeviceVector : public MinlinDeviceVector {
		
		public:
			DeviceVector(int N): MinlinDeviceVector(N) {} /* constructor */
			
	};



}
	



