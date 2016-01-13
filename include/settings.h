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

/* define HostVector/DeviceVector for each data type */
#define DataVector HostVector
#define ThetaVector HostVector
#define GammaVector HostVector
#define GammaMatrix HostMatrix


/* we are going to compute in double/float? */
typedef double Scalar; 
