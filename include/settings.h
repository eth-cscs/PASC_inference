#define DEFAULT_DEBUG_MODE 0

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export results to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */

/* we are going to compute in double/float? */
typedef double Scalar; 

/* define HostVector/DeviceVector for each data type */
#ifdef USE_GPU
	/* compute using CUDA on Device */
/*	#define DataVector DeviceVector<Scalar>
	#define ThetaVector DeviceVector<Scalar>
	#define GammaVector DeviceVector<Scalar>
	#define GammaMatrix DeviceMatrix<Scalar>
*/
#else
	/* compute without CUDA on Host */
/*	#define DataVector HostVector<Scalar>
	#define ThetaVector HostVector<Scalar>
	#define GammaVector HostVector<Scalar>
	#define GammaMatrix HostMatrix<Scalar>
*/
#endif

#ifdef USE_PETSC
	/* compute with Petsc */
	#define DataVector PetscVector
	#define ThetaVector PetscVector
	#define GammaVector PetscVector
	#define GammaMatrix PetscVector


#endif


