#ifndef PETSCVECTOR_H
#define	PETSCVECTOR_H

extern int DEBUG_MODE;
extern bool PETSC_INITIALIZED;

#include <iostream>
#include <string>
#include <list>

static PetscErrorCode ierr;
#define TRY( f) {ierr = f; do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0);}

namespace minlin {

namespace threx { // TODO: maybe choose the different namespace for my own Petsc stuff
 
class PetscVectorWrapperComb; /* wrapper to allow manipulation with linear combinations of vectors */
class PetscVectorWrapperCombNode; /* one node of previous wrapper */
class PetscVectorWrapperSub; /* wrapper to allow subvectors */



/*! \class PetscVector
    \brief Class for manipulation with PetscVector.

    Here will be a detailed description.
*/
class PetscVector {
	private:
		Vec inner_vector; /* original Petsc Vector */
		
	public:

		PetscVector();
		PetscVector(int n);
		PetscVector(const PetscVector &vec1);
		PetscVector(Vec new_inner_vector);
		PetscVector(double scalar_value);
		
		~PetscVector();

		void valuesUpdate();
		void scale(PetscScalar alpha);

		Vec get_vector() const; // TODO: temp, direct access to inner vector should be forbidden
		int size();
		double get(int index);
		void get_ownership(int *low, int *high);

		void get_array(double **arr);
		void restore_array(double **arr);

		void set(double new_value);
		void set(int index, double new_value);

		/* assignment */
		PetscVector &operator=(const PetscVector &vec2);
		PetscVector &operator=(PetscVectorWrapperCombNode combnode);	
		PetscVector &operator=(PetscVectorWrapperComb comb);	
		PetscVector &operator=(double scalar_value);	
		
		/* subvector */
		PetscVectorWrapperSub operator()(int index);
		PetscVectorWrapperSub operator()(int index_begin,int index_end);
		PetscVectorWrapperSub operator()(const IS new_subvector_is);

		friend std::ostream &operator<<(std::ostream &output, const PetscVector &vector);

		friend void operator*=(PetscVector &vec1, double alpha);
		friend void operator+=(PetscVector &vec1, const PetscVectorWrapperComb comb);
		friend void operator-=(PetscVector &vec1, const PetscVectorWrapperComb comb);

		friend double dot(const PetscVector vec1, const PetscVector vec2);
		friend double max(const PetscVector vec1);
		friend double sum(const PetscVector vec1);
		friend const PetscVector operator/(PetscVector vec1, const PetscVector vec2);

		
		/* define operator PetscVector(all), it is a mixture of petscvector and minlin  */
		PetscVectorWrapperSub operator()(minlin::detail::all_type add_in);

};


/*! \class PetscVectorWrapperComb
    \brief Wrapper to allow linear combinations of vectors.

*/
class PetscVectorWrapperComb
{
	private:
		std::list<PetscVectorWrapperCombNode> comb_list; /* the list of linear combination nodes */
		
	public:
		PetscVectorWrapperComb();
		PetscVectorWrapperComb(PetscVectorWrapperCombNode &comb_node);
		PetscVectorWrapperComb(PetscVector &vec);
		PetscVectorWrapperComb(double scalar_value);
		PetscVectorWrapperComb(PetscVectorWrapperSub subvector);

		~PetscVectorWrapperComb();


		int get_listsize();
		int get_vectorsize();
		Vec get_first_vector();
		void append(PetscVectorWrapperCombNode &new_node);
		void merge(PetscVectorWrapperComb &comb);
		void get_arrays(PetscScalar *coeffs, Vec *vectors);
		
		friend std::ostream &operator<<(std::ostream &output, PetscVectorWrapperComb &wrapper);

		friend const PetscVectorWrapperComb operator*(double alpha, PetscVectorWrapperComb comb);
		friend const PetscVectorWrapperComb operator+(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);
		friend const PetscVectorWrapperComb operator-(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);

		
};


/*! \class PetscVectorWrapperComb
    \brief Wrapper with one node in the linear combinations.

*/
class PetscVectorWrapperCombNode
{
	private:
		Vec inner_vector; /* pointer to vector in linear combination */
		double coeff; /* coefficient in linear combination */
	
	public:
		PetscVectorWrapperCombNode(); // TODO: temp
		PetscVectorWrapperCombNode(double new_coeff, Vec new_vector);
		PetscVectorWrapperCombNode(double scalar_value);
	
		void set_vector(Vec new_vector);
		Vec get_vector();
		int get_size();
		int get_value(int index);
		
		void set_coeff(double new_coeff);
		void scale(double alpha);
		double get_coeff();

		friend std::ostream &operator<<(std::ostream &output, PetscVectorWrapperCombNode &wrapper);

};

/*! \class PetscVectorWrapperSub
    \brief Wrapper with subvectors.

*/
class PetscVectorWrapperSub
{
	private:
		IS subvector_is;
		Vec inner_vector; /* original vector */
		Vec subvector; /* subvector */

		bool free_is; /* free index set in destroy */

	public:

		PetscVectorWrapperSub(Vec inner_vector, IS subvector_is, bool new_free_is);
		~PetscVectorWrapperSub();

		void valuesUpdate();
		void scale(PetscScalar alpha);

		void set(double new_value);		
		Vec get_subvector();

		double get(int index);

		friend std::ostream &operator<<(std::ostream &output, const PetscVectorWrapperSub &wrapper);				

		PetscVectorWrapperSub &operator=(const PetscVector &vec2);
		PetscVectorWrapperSub &operator=(PetscVectorWrapperCombNode combnode);	
		PetscVectorWrapperSub &operator=(PetscVectorWrapperComb comb);	
		PetscVectorWrapperSub &operator=(double scalar_value);	

		friend void operator*=(PetscVectorWrapperSub subvec1, double alpha);
		friend void operator+=(PetscVectorWrapperSub subvec1, const PetscVectorWrapperComb comb);
		friend void operator-=(PetscVectorWrapperSub subvec1, const PetscVectorWrapperComb comb);
		friend void operator/=(PetscVectorWrapperSub subvec1, const PetscVectorWrapperSub subvec2);

		friend bool operator==(PetscVectorWrapperSub subvec1, double alpha);
		friend bool operator==(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2);

		friend bool operator>(PetscVectorWrapperSub vec1, PetscVectorWrapperSub vec2);

		friend double sum(const PetscVectorWrapperSub subvec1);
		friend double dot(const PetscVectorWrapperSub subvec1, const PetscVectorWrapperSub subvec2);


};




} /* end of namespace */

} /* end of MinLin namespace */


#include "PetscVector.h"
#include "PetscVectorWrapperComb.h"
#include "PetscVectorWrapperSub.h"

#endif
