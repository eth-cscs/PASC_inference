/** @file petscvector.h
 *  @brief Header file with class delarations.
 *
 *  This defines the basic layout of all used classes. In the end of the file,
 *  the files with specific implementations are included.
 *  This is the file which has to be included in the project to provide the work 
 *  with Petsc vectors in Min-Lin Matlab style.
 *
 *  @author Lukas Pospisil
 *  @bug No known bugs.
 */

#ifndef PETSCVECTOR_H
#define	PETSCVECTOR_H

/* include petsc */
#include "petsc.h"

/* std:list for linear combinations */
#include <list>

/* basic input/output in c++ */
#include <iostream>

/* for manipulating with strings */
#include <string>

/* to deal with errors, call Petsc functions with TRY(fun); */
static PetscErrorCode ierr; /**< to deal with PetscError */

/**
 * \def TRY(f)
 * Macro for dealing with PetscError. Each original Petsc function could by called using TRY(...).
*/
#define TRY( f) {ierr = f; do {if (PetscUnlikely(ierr)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_IN_CXX,0);}} while(0);}

/* we are using namespace petscvector */
namespace petscvector {

int DEBUG_MODE_PETSCVECTOR = true; /**< defines the debug mode of the functions */
bool PETSC_INITIALIZED = false; /**< to deal with PetscInitialize and PetscFinalize outside this class */

/* define "all" stuff */
class petscvector_all_type {} all; /**< brings an opportunity to call PetscVector(all) */

/* wrapper to allow manipulation with linear combinations of vectors */
class PetscVectorWrapperComb;

/* one node of previous wrapper */
class PetscVectorWrapperCombNode;

/* wrapper to allow subvectors */
class PetscVectorWrapperSub; 

/* wrapper to allow (vector or subvector) = mul(v1,v2) */
class PetscVectorWrapperMul; 


/** \class PetscVector
 *  \brief General class for manipulation with vectors.
 *
*/
class PetscVector {
	private:
		Vec inner_vector; /**< original Petsc Vector */
		
	public:

		/** @brief The basic constructor.
		* 
		*  Sets the inner vector to NULL.
		*
		*/
		PetscVector();

		/** @brief Create constructor.
		*
		*  Create new vector of given size n.
		*
		*  @param n global size of new vector
		*/ 
		PetscVector(int n);

		/** @brief Create constructor.
		*
		*  Create sequential vector of given values and size n.
		*
		*  @param values array with values of array
		*  @param n global size of new vector
		*/ 
		PetscVector(double *values, int n);
		
		/** @brief Duplicate constructor.
		*
		*  Create new vector by duplicating given one.
		*
		*  @param vec original vector to be duplicated
		*/ 
		PetscVector(const PetscVector &vec1);

		/** @brief Constructor from Vec.
		*
		*  Construct new vector from given Vec.
		*
		* @param new_inner_vector original Vec
		*/
		PetscVector(const Vec &new_inner_vector);

		/** @brief Constructor from linear combination.
		*
		*  Creates a new vector from given linear combination.
		*
		*  @param new_inner_vector original Vec
		*/ 
		PetscVector(const PetscVectorWrapperComb &comb);

		/** @brief Destructor.
		*
		*  If inner vector is present, then destroy it using VecDestroy.
		*
		*/ 
		~PetscVector();

		/** @brief Update values in inner vector.
		*
		*  Calls VecAssemblyBegin and VecAssemblyEnd.
		*
		*/
		void valuesUpdate() const;

		/** @brief Scale all values in inner vector.
		*
		*  Call VecScale.
		* 
		*  @param alpha the scaling coeficient
		*/ 
		void scale(double alpha);

		/** @brief Get inner vector.
		*
		*  @return original inner vector
		*  @todo this function is temporary
		*/ 
		Vec get_vector() const; // TODO: temp, direct access to inner vector should be forbidden

		/** @brief Get size of inner vector.
		*
		*  @return global size of inner vector
		*  @todo control if inner_vector was allocated
		*/ 
		int size() const;

		/** @brief Get local size of inner vector.
		*
		*  @return local size of inner vector
		*  @todo control if inner_vector was allocated
		*/ 
		int local_size() const;

		/** @brief Get single value.
		*
		*  Return single value with given index of component.
		* 
		*  @note works only with local id, really slow
		*  @return global size of inner vector
		*  @todo control if inner_vector was allocated
		*/ 
		double get(int index);
		
		/** @brief Get ownership of global vector.
		*
		*  Call VecGetOwnershipRange, get the indeces of local components.
		* 
		*  @param low start index
		*  @param high end index + 1
		*  @todo control if inner_vector was allocated
		*/ 
		void get_ownership(int *low, int *high);

		/** @brief Get local array from vector.
		*
		*  Call VecGetArray.
		* 
		*  @note call restore_array after changes in array
		*  @param arr array of vector
		*  @todo control if inner_vector was allocated
		*/ 
 		void get_array(double **arr);

		/** @brief Restore local array to vector.
		*
		*  Call VecRestoreArray.
		* 
		*  @note has to be called after get_array	
		*  @param arr array of vector
		*  @todo control if get_array was called
		*/ 
		void restore_array(double **arr);

		/** @brief Set values in inner vector.
		*
		*  Set all values of the vector to given value, this function is called from overloaded operator.
		*
		*  @param new_value new value of all components
		*  @todo control if inner_vector was allocated
		*/ 
		void set(double new_value);

		/** @brief Update value in inner vector.
		*
		*  Set one specific component of the vector to given value, this function is called from overloaded operator.
		*
		*  @param index index of component
		*  @param new_value new value of component
		*  @todo control if inner_vector was allocated
		*/ 
		void set(int index, double new_value);

		/** @brief Load values from file to PETSC_COMM_SELF.
		*
		*  Uses PetscViewerBinaryOpen, PETSC_COMM_SELF.
		*
		*  @param filename name of file with values
		*  @todo control if file exists
		*/ 
		void load_local(std::string filename);

		/** @brief Load values from file to PETSC_COMM_WORLD.
		*
		*  Uses PetscViewerBinaryOpen, PETSC_COMM_WORLD.
		*
		*  @param filename name of file with values
		*  @todo control if file exists
		*/ 
		void load_global(std::string filename);

		/** @brief Save vector to file in PETSc binary format
		*
		*  Uses PetscViewerBinaryOpen, PETSC_COMM_WORLD.
		*
		*  @param filename name of file
		*  @todo control if vector exists
		*/ 
		void save_binary(std::string filename);

		/** @brief Save vector to file in PETSc ASCII format
		*
		*  Uses PetscViewerASCIIOpen, PETSC_COMM_WORLD.
		*
		*  @param filename name of file
		*  @todo control if vector exists
		*/ 
		void save_ascii(std::string filename);

		/** @brief Assignment operator.
		*
		*  Copy values from one vector to another.
		*  If the inner vector does not exist, then duplicate the vector at first.
		*
		*  @param x vector with new values
		*/ 
		PetscVector &operator=(const PetscVector &x);

		/** @brief Assignment operator.
		*
		*  Set all values in the vector equal to given one.
		*
		*  @param alpha new value
		*  @todo give error if inner_vector does not exist
		*/ 
		PetscVector &operator=(double alpha);

		/** @brief Assignment operator.
		*
		*  Set values of the vector equal to the result from linear combination.
		*  If the inner vector does not exist, then duplicate the vector at first.
		*
		*  @param comb linear combination
		*/ 
		PetscVector &operator=(PetscVectorWrapperComb comb);

		PetscVector &operator=(PetscVectorWrapperMul mul);

		friend void operator*=(PetscVector &vec1, double alpha);
		friend void operator+=(const PetscVector &vec1, const PetscVectorWrapperComb comb);
		friend void operator-=(PetscVector &vec1, const PetscVectorWrapperComb comb);

		/** @brief Get subvector.
		*
		*  Get the subvector with one element defined by given index.
		*
		*  @param index index of the element
		*  @return subvector 
		*/ 
		PetscVectorWrapperSub operator()(int index) const;

		/** @brief Get subvector.
		*
		*  Get the subvector with elements defined by given lower and upper index.
		*
		*  @param index_begin first index of the subvector
		*  @param index_end last index of the subvector
		*  @return subvector 
		*/ 
		PetscVectorWrapperSub operator()(int index_begin,int index_end) const;

		/** @brief Get subvector.
		*
		*  Get the subvector with elements defined by given Petsc Index set.
		*
		*  @param is index set
		*  @return subvector 
		*/ 
		PetscVectorWrapperSub operator()(const IS is) const;

		/** @brief Get subvector.
		*
		*  Get the subvector with all elements.
		*
		*  @param all
		*  @return subvector 
		*/ 
		PetscVectorWrapperSub operator()(petscvector_all_type all) const; /* define operator PetscVector(all)  */

		/** @brief Stream insertion operator.
		*
		*  Prints the content of inner_vector.
		*
		*  @param output output stream
		*  @param vector instance of PetscVector to be printed
		*  @todo print local portion of the vector
		*/ 
		friend std::ostream &operator<<(std::ostream &output, const PetscVector &vector);

		/** @brief Compute dot product.
		*
		*  Computes dot product of two given vectors.
		*  \f[\mathrm{result} = \langle x,y \rangle = \sum\limits_{i = 0}^{size-1} x_i y_i\f]
		*
		*  @param x first vector
		*  @param y second vector
		*  @todo control if inner_vectors were allocated
		*/ 
		friend double dot(const PetscVector &x, const PetscVector &y);

		friend double dot(const PetscVector &x, const PetscVectorWrapperSub y);
		friend double dot(const PetscVectorWrapperSub x, const PetscVector &y);


		/** @brief Get the maximum value in vector.
		*
		*  Computes the maximum value in given vector.
		*  \f[\mathrm{result} = \max \lbrace x_i, i = 0, \dots size-1 \rbrace\f]
		*
		*  @param x vector
		*  @todo control if inner_vector wax allocated
		*/ 
		friend double max(const PetscVector &x);
		
		/** @brief Get the sum of values in vector.
		*
		*  Computes the sum of components of given vector.
		*  \f[\mathrm{result} = \sum\limits_{i = 0}^{size-1} x_i \f]
		*
		*  @param x vector
		*  @todo control if inner_vector wax allocated
		*/ 
		friend double sum(const PetscVector &x);

		/** @brief Get 2-norm of vector.
		*
		*  Computes the sum of components of given vector.
		*  \f[\mathrm{result} = \Vert x \Vert_2 = \sqrt{\sum\limits_{i = 0}^{size-1} x_i^2} \f]
		*
		*  @param x vector
		*  @todo control if inner_vector wax allocated
		*/ 
		friend double norm(const PetscVector &x);

		/** @brief Pointwise divide of two vectors.
		*
		*  Divide values of the inner vector by components of input vector.
		*  \f[ x_i = \frac{x_i}{y_i}, ~~\forall i = 0, \dots, size-1 \f]
		*
		*  @param x vector
		*  @param y vector
		*  @todo this function is really strange 
		*/ 
		friend const PetscVector operator/(const PetscVector &x, const PetscVector &y);

		/** @brief Compute pointwise multiplication.
		*
		*  Computes pointwise product of two given vectors.
		*  \f[\mathrm{mul}_i = x_i y_i \f]
		*  uses Petsc function VecPointwiseMult(Vec w, Vec x,Vec y) 
		* 
		*  @param x first vector
		*  @param y second vector
		*  @todo control if inner_vectors were allocated
		*/ 
		friend PetscVectorWrapperMul mul(const PetscVector &x, const PetscVector &y);

};


/** \class PetscVectorWrapperComb
 *  \brief Wrapper to allow linear combinations of vectors.
 *
 *  The linear combination of vector is stored in the list of PetscVectorWrapperComb.
 *  Finally, using the assignment operator in PetscVector, function compute() is called.
*/
class PetscVectorWrapperComb
{
	private:
		std::list<PetscVectorWrapperCombNode> comb_list; /**< the list of linear combination nodes */
		int vector_size; /**< stores the global size of the last added vector */
		
	public:
		/** @brief The basic constructor.
		* 
		*  @todo deprecated - there is no reason to create empty combination
		*/
		PetscVectorWrapperComb();

		/** @brief Constructor from stand-alone vector.
		* 
		*  If the vector is in given in linear combination without coefficient, 
		*  it will be added to the combination as a node with coefficient 1.0.
		* 
		*  @param vec the vector in linear combination
		*/
		PetscVectorWrapperComb(const PetscVector &vec);

		/** @brief Constructor from combination node.
		* 
		*  Each combination node will be recalled to combination. 
		*  Therefore we can operate (define operators) with only combinations, not with nodes. 
		* 
		*  @param comb_node node of linear combination
		*/
		PetscVectorWrapperComb(const PetscVectorWrapperCombNode &comb_node);

		/** @brief Constructor from stand-alone subvector.
		* 
		*  If the subvector is given in linear combination without coefficient, 
		*  it will be added to the combination as a node with coefficient 1.0.
		* 
		*  @param subvec the subvector in linear combination
		*/
		PetscVectorWrapperComb(PetscVectorWrapperSub subvec);

		/** @brief Destructor.
		* 
		*  Destroy the list with linear combination nodes.
		* 
		*  @todo how to destroy the list? element after element?
		*/
		~PetscVectorWrapperComb();

		/** @brief Get number of nodes in linear combination.
		* 
		*  Get the size of inner list of linear combination nodes.
		* 
		*  @return number of nodes in linear combination
		*/
		int get_listsize() const;

		/** @brief Get the dimension of vectors.
		* 
		*  Get the size of the last added vector.
		* 
		*  @return size of the vector
		*/
		int get_vectorsize() const;

		/** @brief Get the first vector in linear combination.
		* 
		*  Usefull if we want to allocate the vector on the left side of assignment operator,
		*  i.e. if the result vector was not allocated.
		* 
		*  @return vector
		*  @todo private?
		*/
		Vec get_first_vector();
		
		/** @brief Append node to linear combination.
		* 
		*  Append new node to the existing list of linear combination nodes.
		* 
		*  @param node new node to be appended
		*/
		void append(const PetscVectorWrapperCombNode &node);

		/** @brief Merge with given combination.
		* 
		*  Merge two lists with linear combinations.
		*  
		*  @param comb the second combinations
		*/
		void merge(const PetscVectorWrapperComb &comb);

		/** @brief Perform the linear combination.
		* 
		*  Perform the linear combination and store result in given vector.
		*  init_scale = 0 if the method was called from operator=
		*  init_scale = 1 if the method was called from operator+=
		* 
		*  @param y result
		*  @param init_scale initial scale of the result vector
		*/
		void compute(const Vec &y, double init_scale);
		
		/* print */
		friend std::ostream &operator<<(std::ostream &output, PetscVectorWrapperComb comb);

		/** @brief Scale.
		* 
		*  Scale all coefficients of nodes in the linear combination.
		*  
		*  @param alpha scalar
		*  @param comb linear combination to be scaled
		*/
		friend const PetscVectorWrapperComb operator*(double alpha, PetscVectorWrapperComb comb);

		/** @brief Addition operator.
		* 
		*  Add one combination to another.
		*  
		*  @param comb1 first linear combination
		*  @param comb2 second linear combination
		*/
		friend const PetscVectorWrapperComb operator+(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);

		/** @brief Subtraction operator.
		* 
		*  Substract one combination from another, i.e. scale the second one with -1.0 and perform addition.
		*  
		*  @param comb1 first linear combination
		*  @param comb2 second linear combination
		*/
		friend const PetscVectorWrapperComb operator-(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);

		friend const PetscVectorWrapperComb operator+(PetscVectorWrapperComb comb1, double scalar);
		friend const PetscVectorWrapperComb operator+(double scalar,PetscVectorWrapperComb comb2);

};


/*! \class PetscVectorWrapperComb
    \brief Wrapper with one node in the linear combinations.

*/
class PetscVectorWrapperCombNode
{
	private:
		Vec inner_vector; /**< pointer to vector (original Petsc Vec) in linear combination */
		double coeff; /**< coefficient in linear combination */

	public:
		/* constructors and destructor */
		PetscVectorWrapperCombNode();
		PetscVectorWrapperCombNode(const PetscVector &vec);
		PetscVectorWrapperCombNode(double new_coeff, Vec new_vector);
		PetscVectorWrapperCombNode(double new_coeff );

		~PetscVectorWrapperCombNode();


		/* general functions */
		void set_vector(Vec new_vector);
		Vec get_vector() const;
		int get_size() const;
		int get_value(int index) const;
		
		void set_coeff(double new_coeff);
		void scale(double alpha);
		double get_coeff() const;



};

/*! \class PetscVectorWrapperSub
    \brief Wrapper with subvectors.

	For manipulation with subvectors. Based on Petsc function GetSubvector (in constructor) and RestoreSubVector (in destructor).
*/
class PetscVectorWrapperSub
{
	private:
		IS subvector_is; /**< the index set of subvector */
		Vec inner_vector; /**< original vector */
		Vec subvector; /**< subvector (created by VecGetSubVector) */

		bool free_is; /**< free index set in destructor or not */

	public:

		PetscVectorWrapperSub(Vec inner_vector, IS subvector_is, bool new_free_is);
		~PetscVectorWrapperSub();

		void valuesUpdate() const;
		void scale(PetscScalar alpha) const;

		void set(double new_value);		
		Vec get_subvector();

		double get(int index);

		/* print */
		friend std::ostream &operator<<(std::ostream &output, const PetscVectorWrapperSub &wrapper);				

		/* assignment operator */
		PetscVectorWrapperSub &operator=(PetscVectorWrapperSub subvec); /* subvec = subvec */
		PetscVectorWrapperSub &operator=(double scalar_value);	 /* subvec = const */
		PetscVectorWrapperSub &operator=(const PetscVector &vec2); /* subvec = vec */
		PetscVectorWrapperSub &operator=(PetscVectorWrapperComb comb);	

		friend void operator*=(const PetscVectorWrapperSub &subvec1, double alpha); /* subvec = alpha*subvec */
		friend void operator+=(const PetscVectorWrapperSub &subvec1, const PetscVectorWrapperComb comb);
		friend void operator-=(const PetscVectorWrapperSub &subvec1, const PetscVectorWrapperComb comb);
		friend void operator/=(const PetscVectorWrapperSub &subvec1, const PetscVectorWrapperSub subvec2);

		/* boolean operations */
		friend bool operator==(PetscVectorWrapperSub subvec1, double alpha);
		friend bool operator==(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2);

		friend bool operator>(PetscVectorWrapperSub vec1, PetscVectorWrapperSub vec2);

		/* binary operations */
		friend double sum(const PetscVectorWrapperSub subvec1);
		friend double dot(const PetscVectorWrapperSub subvec1, const PetscVectorWrapperSub subvec2);

		friend double dot(const PetscVector &x, const PetscVectorWrapperSub y);
		friend double dot(const PetscVectorWrapperSub x, const PetscVector &y);


		/** @brief Compute pointwise multiplication.
		*
		*  Computes pointwise product of two given subvectors.
		*  \f[\mathrm{mul}_i = x_i y_i \f]
		*  uses Petsc function VecPointwiseMult(Vec w, Vec x,Vec y) 
		* 
		*  @param x first vector
		*  @param y second vector
		*  @todo control if inner_vectors were allocated
		*/ 
		friend PetscVectorWrapperMul mul(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2);

};

/*! \class PetscVectorWrapperMul
    \brief Wrapper for manipulation with mul(v1,v2).
*/
class PetscVectorWrapperMul
{
	private:
		Vec inner_vector1; /**< first vector */
		Vec inner_vector2; /**< second vector */
	public:

		PetscVectorWrapperMul(Vec inner_vector1, Vec innervector2);
		~PetscVectorWrapperMul();

		/** @brief Compute pointwise multiplication.
		*
		*  Computes pointwise product of two given subvectors.
		*  \f[\mathrm{mul}_i = x_i y_i \f]
		*  uses Petsc function VecPointwiseMult(Vec w, Vec x,Vec y) 
		* 
		*  @param result output vector
		*  @todo control if inner_vectors were allocated
		*/ 
		void mul(Vec result);

		Vec get_vector1() const;
		Vec get_vector2() const;
		
};



} /* end of petsc vector namespace */


/* add implementations */
#include "petscvector_impl.h"
#include "wrappercomb_impl.h"
#include "wrappersub_impl.h"
#include "wrappermul_impl.h"

#endif
