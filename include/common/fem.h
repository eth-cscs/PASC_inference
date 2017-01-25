/** @file fem.h
 *  @brief class for reduction and prolongation on fem meshes
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM_H
#define	PASC_FEM_H

/* this class is for petscvector */
typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace common {

/** \class FEM manipulation
 *  \brief Reduction/prolongation between FEM meshes.
 *
*/
class Fem {
	protected:
		Decomposition *decomposition1; /**< decomposition of the larger problem */
		Decomposition *decomposition2; /**< decomposition of smaller problem */
		
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem(Decomposition *decomposition1, Decomposition *decomposition2);

		/** @brief destructor
		*/
		~Fem();
		
		void reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
		void prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

};

/* ----------------- Fem implementation ------------- */

Fem::Fem(Decomposition *decomposition1, Decomposition *decomposition2){
	LOG_FUNC_BEGIN

	this->decomposition1 = decomposition1;
	this->decomposition2 = decomposition2;

	LOG_FUNC_END
}

Fem::~Fem(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

void Fem::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double mysum;
	double diff = (decomposition1->get_T())/(double)(decomposition2->get_T());
	
	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition1->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

		for(int t2 =0; t2 < decomposition2->get_T(); t2++){
			mysum = 0.0;
			for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
				mysum += gammak1_arr[i];
			}
			gammak2_arr[t2] = mysum;
		}

		TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	LOG_FUNC_END
}

void Fem::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	int idx;
	double diff = (decomposition2->get_T()-1.0)/(double)((decomposition1->get_T()-1.0));
	
	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition1->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

		for(int t =0; t < decomposition1->get_T(); t++){
			idx = round(t*diff);
			gammak1_arr[t] = gammak2_arr[idx];
		}

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	TRYCXX( VecRestoreArray(gamma1_Vec,&gammak1_arr) );
	TRYCXX( VecRestoreArray(gamma2_Vec,&gammak2_arr) );

	
	LOG_FUNC_END
}




}
} /* end of namespace */

#endif
