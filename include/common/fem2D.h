/** @file femhat.h
 *  @brief class for reduction and prolongation on fem meshes using hat functions
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM2D_H
#define	PASC_FEM2D_H

#ifndef USE_PETSC
 #error 'FEM2D is for PETSC'
#endif

#include "common/fem.h"
#include "algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace common {

/** \class FEM 2D manipulation
 *  \brief Reduction/prolongation between FEM meshes.
 *
*/
class Fem2D : public Fem {
	protected:
		bool left_overlap;			/**< is there overlap to the left side of space? */
		bool right_overlap;			/**< is there overlap to the right side of space? */
		bool top_overlap;			/**< is there overlap to the top side of space? */
		bool bottom_overlap;		/**< is there overlap to the bottom side of space? */
		
		void compute_overlaps();
		int overlap1_idx_size;		/**< number of elements in overlap1_idx */
		int *overlap1_idx;			/**< permutated indexes of overlap part in grid1 */
		int overlap2_idx_size;		/**< number of elements in overlap2_idx */
		int *overlap2_idx;			/**< permutated indexes of overlap part in grid2 */

		double diff_x;
		double diff_y;

		BGMGraphGrid2D *grid1;
		BGMGraphGrid2D *grid2;

		int *bounding_box1;		/**< bounds of local domain [x1_min,x1_max,y1_min,y1_max] of grid1 */
		int *bounding_box2;		/**< bounds of local domain [x2_min,x2_max,y2_min,y2_max] of grid2 */

	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem2D(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		Fem2D(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem2D();

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		std::string get_name() const;
		
		void reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
		void prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

		void compute_decomposition_reduced();

};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
#endif



/* ----------------- Fem implementation ------------- */
Fem2D::Fem2D(double fem_reduce) : Fem(fem_reduce){
	LOG_FUNC_BEGIN
	
	/* I don't have this information without decompositions */
	this->diff_x = 0;
	this->diff_y = 0;
	
	this->grid1 = NULL;
	this->grid2 = NULL;
	
	this->bounding_box1 = new int[4];
	set_value_array(4, this->bounding_box1, 0); /* initial values */

	this->bounding_box2 = new int[4];
	set_value_array(4, this->bounding_box2, 0); /* initial values */

	
	LOG_FUNC_END
}


Fem2D::Fem2D(Decomposition *decomposition1, Decomposition *decomposition2, double fem_reduce) : Fem(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	this->grid1 = (BGMGraphGrid2D*)(this->decomposition1->get_graph());
	this->grid2 = (BGMGraphGrid2D*)(this->decomposition2->get_graph());

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff = 1; /* time */
	diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);

	if(is_reduced()){
		this->bounding_box1 = new int[4];
		set_value_array(4, this->bounding_box1, 0); /* initial values */
		this->bounding_box2 = new int[4];
		set_value_array(4, this->bounding_box2, 0); /* initial values */

		compute_overlaps();
	}

	LOG_FUNC_END
}

Fem2D::~Fem2D(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

std::string Fem2D::get_name() const {
	return "FEM2D";
}

void Fem2D::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(is_reduced()) << std::endl;
	output_global <<  " - diff             : " << diff << std::endl;
	output_global <<  " - diff_x           : " << diff_x << std::endl;
	output_global <<  " - diff_y           : " << diff_y << std::endl;

	output_global <<  " - bounding_box" << std::endl;
	output_local <<   "   - bounding_box1    : " << print_array(this->bounding_box1, 4) << std::endl;
	output_local <<   "   - bounding_box2    : " << print_array(this->bounding_box2, 4) << std::endl;
	output_local.synchronize();

	output_global <<  " - fem_reduce       : " << fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	if(decomposition1 == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		decomposition1->print(output_global);
		output_global.pop();
	}
	if(grid1 == NULL){
		output_global <<  " - grid1            : NO" << std::endl;
	} else {
		output_global <<  " - grid1            : YES [" << grid1->get_width() << ", " << grid1->get_height() << "]" << std::endl;
		output_global.push();
		grid1->print(output_global);
		output_global.pop();
	}

	if(decomposition2 == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		decomposition2->print(output_global);
		output_global.pop();
	}
	if(grid2 == NULL){
		output_global <<  " - grid2            : NO" << std::endl;
	} else {
		output_global <<  " - grid2            : YES [" << grid2->get_width() << ", " << grid2->get_height() << "]" << std::endl;
		output_global.push();
		grid2->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

void Fem2D::compute_overlaps() {
	LOG_FUNC_BEGIN
	
	if(is_reduced()){
		int width1 = grid1->get_width();
		int height1 = grid1->get_height();
		int width2 = grid2->get_width();
		int height2 = grid2->get_height();
		
		/* get arrays of grids */
		int *DD_affiliation1 = grid1->get_DD_affiliation(); 
		int *DD_permutation1 = grid1->get_DD_permutation(); 
		int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 

		int *DD_affiliation2 = grid2->get_DD_affiliation(); 
		int *DD_permutation2 = grid2->get_DD_permutation(); 
		int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

		/* prepare overlap indexes */
		overlap1_idx_size = (bounding_box1[1]-bounding_box1[0]+1)*(bounding_box1[3]-bounding_box1[2]+1);
		overlap1_idx = new int[overlap1_idx_size];
		overlap2_idx_size = (bounding_box2[1]-bounding_box2[0]+1)*(bounding_box2[3]-bounding_box2[2]+1);
		overlap2_idx = new int[overlap2_idx_size];
		
		/* fill overlapping indexes with.. indexes */
		for(int id_x1 = bounding_box1[0]; id_x1 <= bounding_box1[1]; id_x1++){
			for(int id_y1 = bounding_box1[2]; id_y1 <= bounding_box1[3]; id_y1++){
				overlap1_idx[(id_y1-bounding_box1[2])*(bounding_box1[1]-bounding_box1[0]+1) + (id_x1-bounding_box1[0])] = DD_permutation1[id_y1*width1 + id_x1];
			}
		}
		for(int id_x2 = bounding_box2[0]; id_x2 <= bounding_box2[1]; id_x2++){
			for(int id_y2 = bounding_box2[2]; id_y2 <= bounding_box2[3]; id_y2++){
				overlap2_idx[(id_y2-bounding_box2[2])*(bounding_box2[1]-bounding_box2[0]+1) + (id_x2-bounding_box2[0])] = DD_permutation2[id_y2*width2 + id_x2];
			}
		}
		
	}
	
	LOG_FUNC_END
}

void Fem2D::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;

	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak1_overlap_is;
	Vec gammak1_overlap_Vec;

	int *DD_permutation1 = grid1->get_DD_permutation(); 
	int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 
	int *DD_permutation2 = grid2->get_DD_permutation(); 
	int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

	int Rbegin1 = decomposition1->get_Rbegin();
	int Rbegin2 = decomposition2->get_Rbegin();

	int width1 = grid1->get_width();
	int width2 = grid2->get_width();
	int width_overlap1 = bounding_box1[1] - bounding_box1[0] + 1;
	int height_overlap1 = bounding_box1[3] - bounding_box1[2] + 1;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateGeneral(PETSC_COMM_SELF,overlap1_idx_size,overlap1_idx,PETSC_USE_POINTER,&gammak1_overlap_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int r2=0; r2 < decomposition2->get_Rlocal(); r2++){
				int id2 = DD_permutation2[Rbegin2 + r2];
				int id_y2 = floor(id2/(double)width2);
				int id_x2 = id2 - id_y2*width2;

				/* coordinates in overlap */
				double center_x1 = (id_x2)*diff_x - bounding_box1[0];
				double left_x1 = (id_x2-1)*diff_x - bounding_box1[0];
				double right_x1 = (id_x2+1)*diff_x - bounding_box1[0];

				double center_y1 = (id_y2)*diff_y - bounding_box1[2];
				double left_y1 = (id_y2-1)*diff_y - bounding_box1[2];
				double right_y1 = (id_y2+1)*diff_y - bounding_box1[2];
				
				double mysum = 0.0;
				int counter = 0;
				for(int x1 = floor(left_x1); x1 < right_x1; x1++){
					for(int y1 = floor(left_y1); y1 < right_y1; y1++){
						if(x1 >= 0 && x1 < width_overlap1 && y1 >= 0 && y1 < height_overlap1){
							mysum += gammak1_arr[y1*width_overlap1 + x1];
							counter += 1;
						}
					}
				}
				gammak2_arr[r2] = mysum;// /(double)counter;
				
//				coutAll << "r2 = " << r2 << ", counter = " << counter << ", mysum = " << mysum << std::endl;
			}
//			coutAll.synchronize();

			TRYCXX( VecRestoreArray(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_femhat_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition1->get_Tbegin(), decomposition2->get_Tbegin(), decomposition1->get_Tlocal(), decomposition2->get_Tlocal(), left_t1_idx, left_t2_idx, diff);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );
		TRYCXX( ISDestroy(&gammak1_overlap_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	LOG_FUNC_END
}

void Fem2D::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak2_overlap_is;
	Vec gammak2_overlap_Vec;

	int *DD_permutation1 = grid1->get_DD_permutation(); 
	int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 
	int *DD_permutation2 = grid2->get_DD_permutation(); 
	int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

	int Rbegin1 = decomposition1->get_Rbegin();
	int Rbegin2 = decomposition2->get_Rbegin();

	int width1 = grid1->get_width();
	int height1 = grid1->get_height();
	int width2 = grid2->get_width();
	int height2 = grid2->get_height();
	int width_overlap2 = bounding_box2[1] - bounding_box2[0] + 1;
	int height_overlap2 = bounding_box2[3] - bounding_box2[2] + 1;

	for(int k=0;k<decomposition1->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateGeneral(PETSC_COMM_SELF,overlap2_idx_size,overlap2_idx,PETSC_USE_POINTER,&gammak2_overlap_is) );
		TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_overlap_is, &gammak2_overlap_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_overlap_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int r1=0; r1 < decomposition1->get_Rlocal(); r1++){
				int id1 = DD_invpermutation1[Rbegin1 + r1];
				int id_y1 = floor(id1/(double)width1);
				int id_x1 = id1 - id_y1*width1;

				/* coordinates in overlap */
				int center_x2 = floor((id_x1)/diff_x) - bounding_box2[0];
				int center_y2 = floor((id_y1)/diff_y) - bounding_box2[2];
				
//				gammak1_arr[r1] = GlobalManager.get_rank()/(double)GlobalManager.get_size();
//				gammak1_arr[r1] = id1/((double)(width1*height1));
				gammak1_arr[r1] = gammak2_arr[center_y2*width_overlap2 + center_x2];
				
//				coutAll << "r1 = " << r1 << ", value = " << gammak2_arr[center_y2*width_overlap2 + center_x2] << std::endl;
			}
//			coutAll.synchronize();

			TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_overlap_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_femhat_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, decomposition1->get_T(), decomposition2->get_T(), decomposition1->get_Tbegin(), decomposition2->get_Tbegin(), decomposition1->get_Tlocal(), decomposition2->get_Tlocal(), left_t1_idx, left_t2_idx, diff);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak2_Vec, gammak2_overlap_is, &gammak2_overlap_Vec) );
		TRYCXX( ISDestroy(&gammak2_overlap_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );
	}

	TRYCXX( VecRestoreArray(gamma1_Vec,&gammak1_arr) );
	TRYCXX( VecRestoreArray(gamma2_Vec,&gammak2_arr) );

	LOG_FUNC_END
}

void Fem2D::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	/* decomposition1 has to be set */
	this->grid1 = (BGMGraphGrid2D*)(this->decomposition1->get_graph());

	if(is_reduced()){
		int T_reduced = 1;
		int width_reduced = ceil(grid1->get_width()*fem_reduce);
		int height_reduced = ceil(grid1->get_height()*fem_reduce);
		
		this->grid2 = new BGMGraphGrid2D(width_reduced, height_reduced);
		this->grid2->process_grid();
		
		/* decompose second grid based on the decomposition of the first grid */
		this->grid2->decompose(this->grid1, this->bounding_box1, this->bounding_box2);
		this->grid2->print(coutMaster);
		
		/* compute new decomposition */
		decomposition2 = new Decomposition(T_reduced, 
				*(this->grid2), 
				decomposition1->get_K(), 
				decomposition1->get_xdim(), 
				decomposition1->get_DDT_size(), 
				decomposition1->get_DDR_size());

		compute_overlaps();
	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		set_decomposition_reduced(decomposition1);
		this->grid2 = this->grid1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	diff = 1; /* time */
	diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);
			
	LOG_FUNC_END
}

#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double center_t1 = (Tbegin2+t2)*diff;
		double left_t1 = (Tbegin2+t2-1)*diff;
		double right_t1 = (Tbegin2+t2+1)*diff;
				
		int id_counter = floor(left_t1) - left_t1_idx; /* first index in provided local t1 array */

		double phi_value; /* value of basis function */

		/* left part of hat function */
		double mysum = 0.0;
		int t1 = floor(left_t1);

		/* compute linear combination with coefficients given by basis functions */
		while(t1 <= center_t1){
			phi_value = (t1 - left_t1)/(center_t1 - left_t1);
			mysum += phi_value*data1[id_counter];
			t1 += 1;
			id_counter += 1;
		}

		/* right part of hat function */
		while(t1 < right_t1){
			phi_value = (t1 - right_t1)/(center_t1 - right_t1);
			mysum += phi_value*data1[id_counter];
			t1 += 1;
			id_counter += 1;
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff) {
	int t1 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t1 < T1local){
		int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
		int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

		double t1_left = t2_left_id_orig*diff;
		double t1_right = t2_right_id_orig*diff;

		int t2_left_id = t2_left_id_orig - left_t2_idx;
		int t2_right_id = t2_right_id_orig - left_t2_idx;

		/* value of basis functions */
		double t1_value = 0.0;
		double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
		t1_value += phi_value_left*data2[t2_right_id];
				
		double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
		t1_value += phi_value_right*data2[t2_left_id];

		data1[t1] = t1_value;
	}
}

#endif




}
} /* end of namespace */

#endif
