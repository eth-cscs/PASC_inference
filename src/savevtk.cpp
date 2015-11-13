#include "savevtk.h"

PetscErrorCode save_VTK(Data data, Gamma gamma) 
{
	PetscErrorCode ierr;

	PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	
	
	PetscInt i,k;							/* iterators */
	PetscViewer viewer; 					/* to write to file */
	PetscScalar *x1,*x2; 						/* solution in array, TODO n=3? */

	PetscScalar *gamma_arr,*gamma_max_arr; 						/* gamma in array, maximum numbers of gamma */
	PetscInt *gamma_maxid_arr;				/* indexes of max gamma */
		
	char name_of_file[256];					/* the name of output VTK file */
	
	PetscFunctionBegin;

	/* set MPI variables */
    MPI_Comm_size(PETSC_COMM_WORLD,&proc_n);
    MPI_Comm_rank(PETSC_COMM_WORLD,&proc_id);

	/* write to the name of file the number of processor */
	sprintf(name_of_file, "output/data.vtk");

	/* open file to write */
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, name_of_file, &viewer);CHKERRQ(ierr);
	ierr = PetscViewerASCIISynchronizedAllow(viewer, PETSC_TRUE);CHKERRQ(ierr);

	/* write header to file */
	ierr = PetscViewerASCIIPrintf(viewer,"# vtk DataFile Version 3.1\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"this is my kmeans data with solution\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"ASCII\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"DATASET UNSTRUCTURED_GRID\n");CHKERRQ(ierr);

	/* points - coordinates */
	ierr = PetscViewerASCIIPrintf(viewer,"POINTS %d FLOAT\n", data.get_global_size());CHKERRQ(ierr);
	
	ierr = VecGetArray(data.data_vecs[0],&x1);CHKERRQ(ierr);	
    ierr = VecGetArray(data.data_vecs[1],&x2);CHKERRQ(ierr);	
	for(i=0;i < data.get_local_size();i++){
		ierr = PetscViewerASCIISynchronizedPrintf(viewer,"%f %f %f\n", x1[i], x2[i], 0.0);CHKERRQ(ierr);
	}
	ierr = PetscViewerFlush(viewer); CHKERRQ(ierr);
	ierr = VecRestoreArray(data.data_vecs[0],&x1);CHKERRQ(ierr);	
    ierr = VecRestoreArray(data.data_vecs[1],&x2);CHKERRQ(ierr);	
	ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);

	/* values is points */
	ierr = PetscViewerASCIIPrintf(viewer,"POINT_DATA %d\n", data.get_global_size());CHKERRQ(ierr);

	/* domains */
	ierr = PetscViewerASCIIPrintf(viewer,"SCALARS domain float 1\n",0);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"LOOKUP_TABLE default\n");CHKERRQ(ierr);
	for(i=0;i<data.get_local_size();i++){
		ierr = PetscViewerASCIISynchronizedPrintf(viewer,"%f\n", (PetscScalar)proc_id);CHKERRQ(ierr);
	}
	ierr = PetscViewerFlush(viewer); CHKERRQ(ierr);


	/* prepare gamma data */
	ierr = PetscMalloc(data.get_local_size()*sizeof(PetscScalar),&gamma_max_arr); CHKERRQ(ierr);
	ierr = PetscMalloc(data.get_local_size()*sizeof(PetscInt),&gamma_maxid_arr); CHKERRQ(ierr);
	ierr = VecGetArray(gamma.gamma_vecs[0],&gamma_arr); CHKERRQ(ierr);
	/* write gamma0 */
	ierr = PetscViewerASCIIPrintf(viewer,"SCALARS gamma_%d float 1\n",0);CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"LOOKUP_TABLE default\n");CHKERRQ(ierr);
	for(i=0;i<data.get_local_size();i++){
		gamma_max_arr[i] = gamma_arr[i];
		gamma_maxid_arr[i] = 0;
		ierr = PetscViewerASCIISynchronizedPrintf(viewer,"%f\n", (PetscScalar)gamma_arr[i]);CHKERRQ(ierr);
	}
	ierr = PetscViewerFlush(viewer); CHKERRQ(ierr);
	ierr = VecRestoreArray(gamma.gamma_vecs[0],&gamma_arr); CHKERRQ(ierr);
	/* find the maximum */
	for(k=1;k<gamma.get_dim();k++){
		ierr = VecGetArray(gamma.gamma_vecs[k],&gamma_arr); CHKERRQ(ierr);
		/* write gamma_k */
		ierr = PetscViewerASCIIPrintf(viewer,"SCALARS gamma_%d float 1\n",k);CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(viewer,"LOOKUP_TABLE default\n");CHKERRQ(ierr);
		for(i=0;i<data.get_local_size();i++){
			if(gamma_max_arr[i] < gamma_arr[i]){
				gamma_max_arr[i] = gamma_arr[i];
				gamma_maxid_arr[i] = k;
			}
			ierr = PetscViewerASCIISynchronizedPrintf(viewer,"%f\n", (PetscScalar)gamma_arr[i]);CHKERRQ(ierr);
		}
		ierr = PetscViewerFlush(viewer); CHKERRQ(ierr);
		ierr = VecRestoreArray(gamma.gamma_vecs[k],&gamma_arr); CHKERRQ(ierr);
	}
	/* store the values with max id */
	ierr = PetscViewerASCIIPrintf(viewer,"SCALARS gamma_max_id float 1\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(viewer,"LOOKUP_TABLE default\n");CHKERRQ(ierr);
	for(i=0;i<data.get_local_size();i++){
		ierr = PetscViewerASCIISynchronizedPrintf(viewer,"%f\n", (PetscScalar)gamma_maxid_arr[i]);CHKERRQ(ierr);
	}
	ierr = PetscViewerFlush(viewer); CHKERRQ(ierr);

	/* destroy used stuff */
	ierr = PetscFree(gamma_max_arr); CHKERRQ(ierr);
	ierr = PetscFree(gamma_maxid_arr); CHKERRQ(ierr);

	/* close viewer */
	ierr = PetscViewerASCIISynchronizedAllow(viewer, PETSC_FALSE);CHKERRQ(ierr);	
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);  	

}
