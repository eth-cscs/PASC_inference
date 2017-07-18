#include "external/petscvector/common/decomposition.h"

namespace pascinference {
namespace algebra {

template<>
Decomposition<PetscVector>::Decomposition(int T, int R, int K, int xdim, int DDT_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = R;
	this->K = K;
	this->xdim = xdim;

	this->DDT_size = DDT_size;
	this->DDR_size = 1;

	/* prepare new layout for T */
	destroy_DDT_arrays = true;
	DDT_ranges = new int[this->DDT_size+1];

	Vec DDT_layout;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&DDT_layout) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(DDT_layout, VECMPICUDA));
	#endif

	TRYCXX( VecSetSizes(DDT_layout, PETSC_DECIDE, T ));
	TRYCXX( VecSetFromOptions(DDT_layout) );

	/* get properties of this layout */
	const int *DDT_ranges_const;
	TRYCXX( VecGetOwnershipRanges(DDT_layout, &DDT_ranges_const) );
	for(int i=0;i<DDT_size+1;i++){
		DDT_ranges[i] = DDT_ranges_const[i]; // TODO: how to deal with const int in ranges form PETSc?
	}

	/* destroy temp vector */
	TRYCXX( VecDestroy(&DDT_layout) );

	/* prepare new layout for R */
	/* no graph provided, create one */
	double *coordinates_array = new double[R];
	for(int r=0;r<R;r++){
		coordinates_array[r] = r;
	}
	graph = new BGMGraph<PetscVector>(coordinates_array, R, 1);
	graph->process(0.0);

	destroy_DDR_graph = true; /* I created graph, I should destroy it */
	set_graph(*graph, DDR_size);

	compute_rank();

	LOG_FUNC_END
}

template<>
Decomposition<PetscVector>::Decomposition(int T, BGMGraph<PetscVector> &new_graph, int K, int xdim, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();
	this->K = K;
	this->xdim = xdim;

	/* prepare new layout for R */
	destroy_DDR_graph = false;
	set_graph(new_graph, DDR_size);

	this->DDT_size = 1;
	this->DDR_size = new_graph.get_DD_size();

	/* prepare new layout for T */
	destroy_DDT_arrays = true;
	DDT_ranges = new int[this->DDT_size+1];
	DDT_ranges[0] = 0;
	DDT_ranges[1] = T;

	compute_rank();

	LOG_FUNC_END
}

template<>
Decomposition<PetscVector>::Decomposition(int T, BGMGraph<PetscVector> &new_graph, int K, int xdim, int DDT_size, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();
	this->K = K;
	this->xdim = xdim;

	/* prepare new layout for R */
	destroy_DDR_graph = false;
	set_graph(new_graph, DDR_size);

	this->DDT_size = DDT_size;
	this->DDR_size = graph->get_DD_size();

	/* prepare new layout for T */
	destroy_DDT_arrays = true;
	/* unfortunatelly, we have to compute distribution of T manually */
	int DDT_optimal_local_size = T/(double)DDT_size;
	int DDT_optimal_local_size_residue = T - DDT_optimal_local_size*DDT_size;
	DDT_ranges = new int[this->DDT_size+1];
	DDT_ranges[0] = 0;
	for(int i=0;i<DDT_size;i++){
		DDT_ranges[i+1] = DDT_ranges[i] + DDT_optimal_local_size;
		if(i < DDT_optimal_local_size_residue){
			DDT_ranges[i+1] += 1;
		}
	}

	compute_rank();

	LOG_FUNC_END
}

template<>
Decomposition<PetscVector>::~Decomposition(){
	LOG_FUNC_BEGIN

	if(destroy_DDT_arrays){
		free(DDT_ranges);
	}

	if(destroy_DDR_graph){
		free(graph);
	}

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::compute_rank(){
	LOG_FUNC_BEGIN

	/* get rank of this processor */
	int rank = GlobalManager.get_rank();

	this->DDT_rank = rank/(double)this->DDR_size;
	this->DDR_rank = rank - (this->DDT_rank)*(this->DDR_size);

	/* compute TRbegin */
	if(DDR_size > 1){
        if(DDT_size > 1){
            /* combination of time and space */
            Vec layout;
            TRYCXX( VecCreate(PETSC_COMM_WORLD, &layout) );
            TRYCXX( VecSetSizes(layout,this->get_Tlocal()*this->get_Rlocal(),this->get_T()*this->get_R()) );
            TRYCXX( VecSetFromOptions(layout) );
            TRYCXX( VecGetOwnershipRange( layout, &TRbegin, NULL) );
            TRYCXX( VecDestroy(&layout));
        } else {
            /* only in space */
            TRbegin = get_Rbegin();
        }
	} else {
        /* only in time */
        TRbegin = get_Tbegin();
    }

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::set_graph(BGMGraph<PetscVector> &new_graph, int DDR_size) {
	/* decompose graph */
	this->DDR_size = DDR_size;
	new_graph.decompose(DDR_size);

	this->graph = &new_graph;
	destroy_DDR_graph = false;
	DDR_affiliation = new_graph.get_DD_affiliation();
	DDR_permutation = new_graph.get_DD_permutation();
	DDR_invpermutation = new_graph.get_DD_invpermutation();
	DDR_lengths = new_graph.get_DD_lengths();
	DDR_ranges = new_graph.get_DD_ranges();

}

template<>
void Decomposition<PetscVector>::set_new_graph(BGMGraph<PetscVector> &new_graph, int DDR_size) {

	if(destroy_DDR_graph){
		free(graph);
	}

    set_graph(new_graph,DDR_size);
    destroy_DDR_graph = false;

	compute_rank();
}

template<>
void Decomposition<PetscVector>::createGlobalVec_gamma(Vec *x_Vec) const {
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int K = this->get_K();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );

	#ifdef USE_CUDA
		TRYCXX(VecSetType(*x_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(*x_Vec, VECMPI));
	#endif

	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*K,T*R*K) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::createGlobalVec_data(Vec *x_Vec) const {
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int xdim = this->get_xdim();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );

	#ifdef USE_CUDA
		TRYCXX(VecSetType(*x_Vec,VECMPICUDA));
	#else
		TRYCXX(VecSetType(*x_Vec,VECMPI));
	#endif

	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*xdim,T*R*xdim) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::permute_TRb_to_dTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const {
	LOG_FUNC_BEGIN

	int TRbegin = get_TRbegin();
	int Tlocal = get_Tlocal();
	int Rlocal = get_Rlocal();
	int local_size = Tlocal*Rlocal*blocksize;

	IS orig_local_is;
	IS new_local_is;

	Vec orig_local_Vec;
	Vec new_local_Vec;

	/* prepare index set with local data */
	int *orig_local_arr;
	orig_local_arr = new int [local_size];

	/* original data is TRb, therefore there is not necessity to transfer it to TRb */
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, local_size, TRbegin, 1 , &orig_local_is) );

	createIS_dTRb(&new_local_is, blocksize);

	/* get subvector with local values from original data */
	TRYCXX( VecGetSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecGetSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* copy values */
	if(!invert){
		TRYCXX( VecCopy(orig_local_Vec, new_local_Vec) );
	} else {
		TRYCXX( VecCopy(new_local_Vec, orig_local_Vec) );
	}

	/* restore subvector with local values from original data */
	TRYCXX( VecRestoreSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecRestoreSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* destroy used stuff */
	TRYCXX( ISDestroy(&orig_local_is) );
	TRYCXX( ISDestroy(&new_local_is) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::permute_TbR_to_dTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const {
	LOG_FUNC_BEGIN

	int Tbegin = get_Tbegin();
	int Tend = get_Tend();
	int Tlocal = get_Tlocal();
	int Rbegin = get_Rbegin();
	int Rend = get_Rend();
	int Rlocal = get_Rlocal();
	int local_size = Tlocal*Rlocal*blocksize;

	IS orig_local_is;
	IS new_local_is;

	Vec orig_local_Vec;
	Vec new_local_Vec;

	/* prepare index set with local data */
	int *orig_local_arr;
	orig_local_arr = new int [local_size];


	/* original data is bTR */
	/* transfer it to TRb */
	for(int t=0;t<Tlocal;t++){
		for(int r=0;r<Rlocal;r++){
            for(int k=0;k<blocksize;k++){
				orig_local_arr[t*Rlocal*blocksize + r*blocksize + k] = (Tbegin+t)*blocksize*R + k*R + (Rbegin+r);
			}
		}
	}

	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_COPY_VALUES,&orig_local_is) );
//	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_OWN_POINTER,&orig_local_is) );

	createIS_dTRb(&new_local_is, blocksize);

	/* get subvector with local values from original data */
	TRYCXX( VecGetSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecGetSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* copy values */
	if(!invert){
		TRYCXX( VecCopy(orig_local_Vec, new_local_Vec) );
	} else {
		TRYCXX( VecCopy(new_local_Vec, orig_local_Vec) );
	}

	/* restore subvector with local values from original data */
	TRYCXX( VecRestoreSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecRestoreSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* destroy used stuff */
	TRYCXX( ISDestroy(&orig_local_is) );
	TRYCXX( ISDestroy(&new_local_is) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::permute_bTR_to_dTRb(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const {
	LOG_FUNC_BEGIN

	int Tbegin = get_Tbegin();
	int Tend = get_Tend();
	int Tlocal = get_Tlocal();
	int Rbegin = get_Rbegin();
	int Rend = get_Rend();
	int Rlocal = get_Rlocal();
	int local_size = Tlocal*Rlocal*blocksize;

	IS orig_local_is;
	IS new_local_is;

	Vec orig_local_Vec;
	Vec new_local_Vec;

	/* prepare index set with local data */
	int *orig_local_arr;
	orig_local_arr = new int [local_size];

	/* original data is TbR */
	/* transfer it to TRb */
	for(int t=0;t<Tlocal;t++){
		for(int r=0;r<Rlocal;r++){
            for(int k=0;k<blocksize;k++){
				orig_local_arr[t*Rlocal*blocksize + r*blocksize + k] = k*R*T + (Tbegin+t)*R + (Rbegin+r);
			}
		}
	}

	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_COPY_VALUES,&orig_local_is) );
//	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_OWN_POINTER,&orig_local_is) );

	createIS_dTRb(&new_local_is, blocksize);

	/* get subvector with local values from original data */
	TRYCXX( VecGetSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecGetSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* copy values */
	if(!invert){
		TRYCXX( VecCopy(orig_local_Vec, new_local_Vec) );
	} else {
		TRYCXX( VecCopy(new_local_Vec, orig_local_Vec) );
	}

	/* restore subvector with local values from original data */
	TRYCXX( VecRestoreSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecRestoreSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* destroy used stuff */
	TRYCXX( ISDestroy(&orig_local_is) );
	TRYCXX( ISDestroy(&new_local_is) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::createIS_dTRb(IS *is, int blocksize) const {
	/* from dTRb to TRb */

	LOG_FUNC_BEGIN

	int Tbegin = get_Tbegin();
	int Tend = get_Tend();
	int Tlocal = get_Tlocal();
	int Rbegin = get_Rbegin();
	int Rend = get_Rend();
	int Rlocal = get_Rlocal();

	int local_size = Tlocal*Rlocal*blocksize;

	/* fill orig_local_arr */
	int *local_arr;
	local_arr = new int[local_size];
	int *DD_permutation = get_DDR_permutation();
	for(int t=0;t<Tlocal;t++){
		for(int r=0;r<Rlocal;r++){
			for(int k=0;k<blocksize;k++){
                local_arr[t*Rlocal*blocksize + r*blocksize + k] = (Tbegin+t)*R*blocksize + DD_permutation[Rbegin+r]*blocksize + k;
            }
		}
	}

	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, local_arr, PETSC_COPY_VALUES,is) );

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::createIS_gammaK(IS *is, int k) const {
	LOG_FUNC_BEGIN

	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, get_Tlocal()*get_Rlocal(), get_TRbegin()*get_K() + k, get_K(), is) );

	LOG_FUNC_END
}

template<>
void Decomposition<PetscVector>::createIS_datan(IS *is, int n) const {
	LOG_FUNC_BEGIN

	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, get_Tlocal()*get_Rlocal(), get_TRbegin()*get_xdim() + n, get_xdim(), is) );

	LOG_FUNC_END
}


}
} /* end of namespace */

