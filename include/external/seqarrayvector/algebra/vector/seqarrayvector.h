/** @file seqarrayvector.h
 *  @brief Header file with class delarations.
 *
 *  @author Lukas Pospisil
 */

#ifndef SEQARRAYVECTOR_H
#define	SEQARRAYVECTOR_H

/* std:list for linear combinations */
#include <list>

/* basic input/output in c++ */
#include <iostream>

/* for manipulating with strings */
#include <string>

#include <limits>
#include <math.h>

/* we are using namespace petscvector */
namespace seqarrayvector {

/* define "all" stuff */
class seqarrayvector_all_type {};
extern seqarrayvector_all_type all; /**< brings an opportunity to call SeqArrayVector(all) */

/** \class SeqArrayVector
 *  \brief General class for manipulation with vectors.
 *
*/
class SeqArrayVector {
	private:
		double *inner_array; /**< content of values */
		int inner_size;		 /**< length of array */
		
	public:

		/** @brief The basic constructor.
		* 
		*  Sets the inner array to NULL.
		*
		*/
		SeqArrayVector();

		/** @brief Create constructor.
		*
		*  Create new vector of given size n.
		*
		*  @param n global size of new vector
		*/ 
		SeqArrayVector(int n);

		/** @brief Create constructor.
		*
		*  Create sequential vector of given values and size n.
		*
		*  @param values array with values of array
		*  @param n global size of new vector
		*/ 
		SeqArrayVector(double *values, int n);
		
		/** @brief Duplicate constructor.
		*
		*  Create new vector by duplicating given one.
		*
		*  @param vec original vector to be duplicated
		*/ 
		SeqArrayVector(const SeqArrayVector &vec1);

		/** @brief Destructor.
		*
		*  If inner array is present, then destroy it.
		*
		*/ 
		~SeqArrayVector();

		/** @brief Scale all values in inner vector.
		*
		*  @param alpha the scaling coeficient
		*/ 
		void scale(double alpha);

		/** @brief Get inner array.
		*
		*  @return original inner array
		*/ 
		double *get_array() const;

		/** @brief Get size of inner array.
		*
		*  @return global size of inner array
		*/ 
		int size() const;

		/** @brief Get local size of inner array.
		*
		*  @return local size of inner vector
		*/ 
		int local_size() const;

		/** @brief Get single value.
		*
		*  Return single value with given index of component.
		* 
		*  @return component of vector
		*/ 
		double get(int index);
		
		/** @brief Get ownership of global vector.
		*
		*  @param low start index
		*  @param high end index + 1
		*/ 
		void get_ownership(int *low, int *high);

		/** @brief Set values in inner vector.
		*
		*  Set all values of the vector to given value, this function is called from overloaded operator.
		*
		*  @param new_value new value of all components
		*/ 
		void set(double new_value);

		/** @brief Update value in inner vector.
		*
		*  Set one specific component of the vector to given value, this function is called from overloaded operator.
		*
		*  @param index index of component
		*  @param new_value new value of component
		*/ 
		void set(int index, double new_value);

		/** @brief Load values from CSV file.
		*
		*  @param filename name of file with values
		*  @todo control if file exists
		*/ 
		void load_csv(std::string filename);

		/** @brief Save vector to file in CSV format
		*
		*  @param filename name of file
		*  @todo control if vector exists
		*/ 
		void save_csv(std::string filename);

		/** @brief Assignment operator.
		*
		*  Copy values from one vector to another.
		*  If the inner vector does not exist, then duplicate the vector at first.
		*
		*  @param x vector with new values
		*/ 
		SeqArrayVector &operator=(const SeqArrayVector &x);

		/** @brief Assignment operator.
		*
		*  Set all values in the vector equal to given one.
		*
		*  @param alpha new value
		*  @todo give error if inner_vector does not exist
		*/ 
		SeqArrayVector &operator=(double alpha);

		friend void operator*=(SeqArrayVector &vec1, double alpha);
		friend void operator+=(const SeqArrayVector &vec1, const SeqArrayVector comb);
		friend void operator-=(SeqArrayVector &vec1, const SeqArrayVector comb);

		/** @brief Get subvector.
		*
		*  Get the subvector with one element defined by given index.
		*
		*  @param index index of the element
		*  @return subvector 
		*/ 
		SeqArrayVector operator()(int index) const;

		/** @brief Get subvector.
		*
		*  Get the subvector with elements defined by given lower and upper index.
		*
		*  @param index_begin first index of the subvector
		*  @param index_end last index of the subvector
		*  @return subvector 
		*/ 
		SeqArrayVector operator()(int index_begin,int index_end) const;

		/** @brief Get subvector.
		*
		*  Get the subvector with all elements.
		*
		*  @param all
		*  @return subvector 
		*/ 
		SeqArrayVector operator()(seqarrayvector_all_type all) const;

		/** @brief Stream insertion operator.
		*
		*  Prints the content of inner array.
		*
		*  @param output output stream
		*  @param vector instance of SeqArrayVector to be printed
		*/ 
		friend std::ostream &operator<<(std::ostream &output, const SeqArrayVector &vector);
		
		/** @brief Compute dot product.
		*
		*  Computes dot product of two given vectors.
		*  \f[\mathrm{result} = \langle x,y \rangle = \sum\limits_{i = 0}^{size-1} x_i y_i\f]
		*
		*  @param x first vector
		*  @param y second vector
		*/ 
		friend double dot(const SeqArrayVector &x, const SeqArrayVector &y);

		/** @brief Get the maximum value in vector.
		*
		*  Computes the maximum value in given vector.
		*  \f[\mathrm{result} = \max \lbrace x_i, i = 0, \dots size-1 \rbrace\f]
		*
		*  @param x vector
		*/ 
		friend double max(const SeqArrayVector &x);

		/** @brief Get the minimum value in vector.
		*
		*  Computes the minimum value in given vector.
		*  \f[\mathrm{result} = \min \lbrace x_i, i = 0, \dots size-1 \rbrace\f]
		*
		*  @param x vector
		*/ 
		friend double min(const SeqArrayVector &x);
		
		/** @brief Get the sum of values in vector.
		*
		*  Computes the sum of components of given vector.
		*  \f[\mathrm{result} = \sum\limits_{i = 0}^{size-1} x_i \f]
		*
		*  @param x vector
		*/ 
		friend double sum(const SeqArrayVector &x);

		/** @brief Get 2-norm of vector.
		*
		*  Computes the sum of components of given vector.
		*  \f[\mathrm{result} = \Vert x \Vert_2 = \sqrt{\sum\limits_{i = 0}^{size-1} x_i^2} \f]
		*
		*  @param x vector
		*/ 
		friend double norm(const SeqArrayVector &x);

		/** @brief Pointwise divide of two vectors.
		*
		*  Divide values of the inner vector by components of input vector.
		*  \f[ x_i = \frac{x_i}{y_i}, ~~\forall i = 0, \dots, size-1 \f]
		*
		*  @param x vector
		*  @param y vector
		*  @todo this function is really strange 
		*/ 
		friend const SeqArrayVector operator/(const SeqArrayVector &x, const SeqArrayVector &y);

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
		friend SeqArrayVector mul(const SeqArrayVector &x, const SeqArrayVector &y);

};

extern void operator*=(SeqArrayVector &vec1, double alpha);
extern std::ostream &operator<<(std::ostream &output, const SeqArrayVector &vector);
extern void operator*=(SeqArrayVector &vec1, double alpha);
extern void operator+=(const SeqArrayVector &vec1, const SeqArrayVector vec2);
extern void operator-=(const SeqArrayVector &vec1, const SeqArrayVector vec2);
extern std::ostream &operator<<(std::ostream &output, const SeqArrayVector &vector);
extern double dot(const SeqArrayVector &x, const SeqArrayVector &y);
extern double max(const SeqArrayVector &x);
extern double min(const SeqArrayVector &x);
extern double sum(const SeqArrayVector &x);
extern double norm(const SeqArrayVector &x);
extern const SeqArrayVector operator/(const SeqArrayVector &x, const SeqArrayVector &y);
extern SeqArrayVector mul(const SeqArrayVector &x, const SeqArrayVector &y);


} /* end of petsc vector namespace */


#endif
