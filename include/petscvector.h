#ifndef PETSCVECTOR_H
#define	PETSCVECTOR_H

#include <iostream>
#include <string>

class PetscVectorWrapperAssign; /* wrapper to allow vector(i) = value */

/* PETSc Vector */
class PetscVector {
		PetscErrorCode ierr; // TODO: I don't know what to do with errors
		Vec inner_vector;
	public:
		/* constructor with given dimension */
		PetscVector(int n){
			ierr = VecCreate(PETSC_COMM_WORLD,&inner_vector);
			ierr = VecSetSizes(inner_vector,PETSC_DECIDE,n);
			ierr = VecSetFromOptions(inner_vector);

			valuesUpdate();
		}

		/* constructor with given internal Vec */
		PetscVector(Vec new_inner_vector){
			inner_vector = new_inner_vector;
		}


		/* destructor */
		~PetscVector(){
//			ierr = VecDestroy(&vector);
		}
		
		/* get PETSC original vector */
		Vec get_vector() const { // TODO: temp
			return inner_vector;
		}
		
		/* get size of the vector */
		int get_size() const {
			int global_size;
			VecGetSize(inner_vector,&global_size);
			return global_size;
		}

		/* after update a variable, it is necessary to call asseble begin */
		void valuesUpdate(){
			VecAssemblyBegin(inner_vector);
			VecAssemblyEnd(inner_vector);
		}
		
		/* set value of the vector, this function is called from overloaded operator */
		void set(int index, PetscScalar new_value){
			VecSetValue(inner_vector,index,new_value, INSERT_VALUES);
			valuesUpdate();
		}

		/* set all values of the vector, this function is called from overloaded operator */
		void set(PetscScalar new_value){
			VecSet(inner_vector,new_value);
			valuesUpdate();
		}

		/* inner_vector = alpha*inner_vector */
		void scale(PetscScalar alpha){
			VecScale(inner_vector, alpha);
			valuesUpdate(); // TODO: has to be called?
		}


		/* stream insertion operator */
		friend std::ostream &operator<<(std::ostream &output, const PetscVector &vector)		
		{
			PetscScalar *arr_vector;
			PetscInt i,local_size;
	
			VecGetLocalSize(vector.inner_vector,&local_size);
			VecGetArray(vector.inner_vector,&arr_vector);
			for (i=0; i<local_size; i++){
				output << arr_vector[i];
				if(i < local_size-1) output << ", ";
			}
			VecRestoreArray(vector.inner_vector,&arr_vector);
			
			return output;
		}

		/* get value with given id of the vector (works only with local id) */
		PetscScalar operator ()(int i) const
		{
			PetscInt ni = 1;
			PetscInt ix[1];
			PetscScalar y[1];
			
			ix[0] = i;
			VecGetValues(inner_vector,ni,ix,y);			
			
			return y[0];
		}

		/* set value with given id of the vector (works only with local id), will be defined after PetscVector */
		PetscVectorWrapperAssign operator()(int index);

		/* vec1 *= alpha */
		void operator*=(double alpha)
		{
			scale(alpha);
		}

		/* assignment operator (copy) */
		PetscVector operator=(const PetscVector &new_vector)
		{
			VecCopy(new_vector.inner_vector,inner_vector);
		}
	
		/* vec1 = alpha*vec2 */
		friend const PetscVector operator*(double alpha, const PetscVector vec2);

		/* vec1 += vec2 */
		friend const void operator+=(PetscVector vec1, const PetscVector vec2);

		/* vec1 -= vec2 */
		friend const void operator-=(PetscVector vec1, const PetscVector vec2);

		/* vec1 = vec2 + vec3 */
		friend const PetscVector operator+(const PetscVector vec2, const PetscVector vec3);

		/* vec1 = vec2 - vec3 */
		friend const PetscVector operator-(const PetscVector vec2, const PetscVector vec3);

		/* dot = dot(vec1,vec2) */
		friend const double dot(const PetscVector vec1, const PetscVector vec2);
	
};

/* wrapper to allow vector(i) = value */
class PetscVectorWrapperAssign
{
	PetscVector store; /* in this vector we want to store new value */
	int index; /* index of new value */
	
	public:
		/* constructor */
		PetscVectorWrapperAssign(PetscVector &s, int i): store(s), index(i) {}
		
		/* define assigment operator */
		PetscVectorWrapperAssign& operator=(PetscScalar const& new_value)
		{
			/* I am not able to access private vector, I pass it to orig class */
			store.set(index,new_value);
		}
	
};

/* return wrapper to be able to overload vector(index) = new_value */ 
PetscVectorWrapperAssign PetscVector::operator()(int index)
{   
	return PetscVectorWrapperAssign(*this, index);
}

/* vec1 = alpha*vec2 */
const PetscVector operator*(double alpha, const PetscVector vec2) // TODO: make a wrapper for linear combinations
{
	Vec new_inner_vector;
	VecDuplicate(vec2.inner_vector,&new_inner_vector);
	VecCopy(vec2.inner_vector,new_inner_vector);
	VecScale(new_inner_vector,alpha);

	return PetscVector(new_inner_vector);
}

/* vec1 += vec2 */
const void operator+=(PetscVector vec1, const PetscVector vec2)
{
	VecAXPY(vec1.inner_vector,1.0, vec2.inner_vector);
}

/* vec1 -= vec2 */
const void operator-=(PetscVector vec1, const PetscVector vec2)
{
	VecAXPY(vec1.inner_vector,-1.0, vec2.inner_vector);
}

/* vec1 = vec2 + vec3 */
const PetscVector operator+(const PetscVector vec2, const PetscVector vec3) // TODO: make a wrapper for linear combinations
{
	Vec new_inner_vector;
	VecDuplicate(vec2.inner_vector,&new_inner_vector); 
	VecCopy(vec2.inner_vector,new_inner_vector); /* vec1 = vec2 */
	VecAXPY(new_inner_vector,1.0, vec3.inner_vector); /* vec1 += vec3 */

	return PetscVector(new_inner_vector);
}

/* vec1 = vec2 - vec3 */
const PetscVector operator-(const PetscVector vec2, const PetscVector vec3) // TODO: make a wrapper for linear combinations
{
	Vec new_inner_vector;
	VecDuplicate(vec2.inner_vector,&new_inner_vector); 
	VecCopy(vec2.inner_vector,new_inner_vector); /* vec1 = vec2 */
	VecAXPY(new_inner_vector,-1.0, vec3.inner_vector); /* vec1 -= vec3 */

	return PetscVector(new_inner_vector);
}

/* dot = dot(vec1,vec2) */
const double dot(const PetscVector vec1, const PetscVector vec2)
{
	double dot_value;
	VecDot(vec1.inner_vector,vec2.inner_vector,&dot_value);
	return dot_value;
}


#endif
