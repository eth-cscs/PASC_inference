/* dimension of the data */
#define datan 2 // number of data components, my sample is in 2D, TODO: do it in the different way

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export results to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */

/* define HostVector/DeviceVector for each data type */
#define DataVector DeviceVector
#define ThetaVector DeviceVector
#define GammaVector DeviceVector
#define GammaMatrix DeviceMatrix


/* we are going to compute in double/float? */
typedef double Scalar; 
