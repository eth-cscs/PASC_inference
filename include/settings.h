/* dimension of the data */
#define datan 2 // number of data components, my sample is in 2D, TODO: do it in the different way

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export results to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */

#define DEBUG_PRINTDATA true /* print values of all data */
#define DEBUG_PRINTL true /* print descent of object function in main outer loop */

#define DEFAULT_T 10 /* default length of generated time serie */
#define DEFAULT_K 3 /* default number of clusters */

#define ALGORITHM_deltaL_eps 0.0001 /*stopping criteria of outer main loop */
#define ALGORITHM_max_s_steps 1 /* max number of outer steps */

#define ALGORITHM_SPGQP_m 10
#define ALGORITHM_SPGQP_gamma 0.9
#define ALGORITHM_SPGQP_sigma2 0.9999
#define ALGORITHM_SPGQP_eps 0.0001
#define ALGORITHM_SPGQP_maxit 10
#define ALGORITHM_SPGQP_alphainit 0.25
#define DEBUG_ALGORITHM_PRINTFS true /* print vector of object functions in every iteration */
#define DEBUG_ALGORITHM_PRINTCOEFF true /* print computed coefficients in every iteration */

/* define HostVector/DeviceVector for each data type */
#define DataVector HostVector
#define ThetaVector HostVector
#define GammaVector HostVector
#define GammaMatrix HostMatrix


/* we are going to compute in double/float? */
typedef double Scalar; 
