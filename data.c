/*
 * data.c
 *
 *  Created on: 18 Jan 2016
 *      Author: ofjp
 */

#include "data.h"

void init_data(struct data_obj *data){

	(*data).x_min = NULL;
	(*data).x_max = NULL;
	(*data).y_min = NULL;
	(*data).y_max = NULL;

	(*data).dt = NULL;
	(*data).cellx = NULL;
	(*data).celly = NULL;
	(*data).vertexdx = NULL;
	(*data).vertexdy = NULL;
	(*data).celldx = NULL;
	(*data).celldy = NULL;
	(*data).xarea = NULL;
	(*data).yarea = NULL;
	(*data).volume = NULL;
	(*data).density0 = NULL;
	(*data).density1 = NULL;
	(*data).energy0 = NULL;
	(*data).energy1 = NULL;
	(*data).pressure = NULL;
	(*data).viscosity = NULL;
	(*data).soundspeed = NULL;
	(*data).xvel0 = NULL;
	(*data).yvel0 = NULL;
	(*data).xvel1 = NULL;
	(*data).yvel1 = NULL;
	(*data).vol_flux_x = NULL;
	(*data).vol_flux_y = NULL;
	(*data).mass_flux_x = NULL;
	(*data).mass_flux_y = NULL;
	(*data).work_array1 = NULL;
	(*data).work_array2 = NULL;
	(*data).work_array3 = NULL;
	(*data).work_array4 = NULL;
	(*data).work_array5 = NULL;
	(*data).work_array6 = NULL;
	(*data).work_array7 = NULL;

	(*data).vertexx = NULL;
	(*data).vertexy = NULL;


}

void allocate_data(struct data_obj data){

    int x_min, x_max, y_min, y_max;
    int j,k;

	printf("Allocating Data\n");

	if((data.x_min == NULL) || (data.x_max == NULL) ||(data.y_min == NULL) ||(data.y_max == NULL) ){
		printf("Error: Allocation with no problem size specified\n");
        return;
	}


	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	int size_xsmall_2d = sizeof(double)*(((x_max-x_min)+4+1)); //(x_min-2:x_max+2)
	int size_xlarge_2d = sizeof(double)*(((x_max-x_min)+5+1)); //(x_min-2:x_max+3)
	int size_ysmall_2d = sizeof(double)*(((y_max-y_min)+4+1)); //(y_min-2:y_max+2)
	int size_ylarge_2d = sizeof(double)*(((y_max-y_min)+5+1)); //(y_min-2:y_max+3)


	int size_small_3d = sizeof(double)*(((x_max-x_min)+4+1)*((y_max-y_min)+4+1)); //(x_min-2:x_max+2,y_min-2:y_max+2)
	int size_x_3d = sizeof(double)*(((x_max-x_min)+5+1)*((y_max-y_min)+4+1)); //(x_min-2:x_max+3,y_min-2:y_max+2)
	int size_y_3d = sizeof(double)*(((x_max-x_min)+4+1)*((y_max-y_min)+5+1)); //(x_min-2:x_max+2,y_min-2:y_max+3)
	int size_large_3d = sizeof(double)*(((x_max-x_min)+5+1)*((y_max-y_min)+5+1)); //(x_min-2:x_max+3,y_min-2:y_max+3)


	if(data.vertexx != NULL){
	    *(data.vertexx)=(double *) malloc(size_xlarge_2d);
	}
	if(data.vertexy != NULL){
		*(data.vertexy)=(double *) malloc(size_ylarge_2d);
	}


	if(data.cellx != NULL){
		*(data.cellx) = (double *) malloc(size_xsmall_2d);
	}

	if(data.celly != NULL){
		*(data.celly) = (double *) malloc(size_ysmall_2d);
	}

	if(data.vertexdx != NULL){
		*(data.vertexdx) = (double *) malloc(size_xlarge_2d);
	}

	if(data.vertexdy != NULL){
		*(data.vertexdy) = (double *) malloc(size_ylarge_2d);
	}


	if(data.celldx != NULL){
		*(data.celldx) = (double *) malloc(size_xsmall_2d);
	}

	if(data.celldy != NULL){
		*(data.celldy) = (double *) malloc(size_ysmall_2d);
	}


	if(data.xarea != NULL){
		*(data.xarea) = (double *) malloc(size_x_3d);
        }

	if(data.yarea != NULL){
		*(data.yarea) = (double *) malloc(size_y_3d);
	}


	if(data.volume != NULL){
		*(data.volume) = (double *) malloc(size_small_3d);
	}


	if(data.density0 != NULL){
		*(data.density0) = (double *) malloc(size_small_3d);
	}
	if(data.density1 != NULL){
		*(data.density1) = (double *) malloc(size_small_3d);
        }

	if(data.energy0 != NULL){
		*(data.energy0) = (double *) malloc(size_small_3d);
        }

	if(data.energy1 != NULL){
		*(data.energy1) = (double *) malloc(size_small_3d);
        }

	if(data.pressure != NULL){
		*(data.pressure) = (double *) malloc(size_small_3d);
        }

	if(data.soundspeed != NULL){
		*(data.soundspeed) = (double *) malloc(size_small_3d);
	}

	if(data.viscosity != NULL){
		*(data.viscosity) = (double *) malloc(size_small_3d);
	}

	if(data.xvel0 != NULL){
		*(data.xvel0) = (double *) malloc(size_large_3d);
	}

	if(data.xvel1 != NULL){
		*(data.xvel1) = (double *) malloc(size_large_3d);
	}

	if(data.yvel0 != NULL){
		*(data.yvel0) = (double *) malloc(size_large_3d);
	}

	if(data.yvel1 != NULL){
		*(data.yvel1) = (double *) malloc(size_large_3d);
	}



	if(data.work_array1 != NULL){
		*(data.work_array1) = (double *) malloc(size_large_3d);
	}
	if(data.work_array2 != NULL){
		*(data.work_array2) = (double *) malloc(size_large_3d);
        }
	if(data.work_array3 != NULL){
		*(data.work_array3) = (double *) malloc(size_large_3d);
	}
	if(data.work_array4 != NULL){
		*(data.work_array4) = (double *) malloc(size_large_3d);
	}
	if(data.work_array5 != NULL){
		*(data.work_array5) = (double *) malloc(size_large_3d);
	}
	if(data.work_array6 != NULL){
		*(data.work_array6) = (double *) malloc(size_large_3d);
	}
	if(data.work_array7 != NULL){
		*(data.work_array7) = (double *) malloc(size_large_3d);
	}


	if(data.vol_flux_x != NULL){
		*(data.vol_flux_x) = (double *) malloc(size_x_3d);
	}
	if(data.vol_flux_y != NULL){
		*(data.vol_flux_y) = (double *) malloc(size_y_3d);
	}

	if(data.mass_flux_x != NULL){
		*(data.mass_flux_x) = (double *) malloc(size_x_3d);
	}
	if(data.mass_flux_y != NULL){
		*(data.mass_flux_y) = (double *) malloc(size_y_3d);
	}
	
#pragma omp parallel
{

	if(data.xarea != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.xarea))[index] = 0.0;
			}
		}
        }

	if(data.yarea != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.yarea))[index] = 0.0;
			}
		}
	}


	if(data.volume != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.volume))[index] = 0.0;
			}
		}
	}


	if(data.density0 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.density0))[index] = 0.0;
			}
		}
	}
	if(data.density1 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.density1))[index] = 0.0;
			}
		}
        }

	if(data.energy0 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.energy0))[index] = 0.0;
			}
		}
        }

	if(data.energy1 != NULL){
		*(data.energy1) = (double *) malloc(size_small_3d);
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.energy1))[index] = 0.0;
			}
		}
        }

	if(data.pressure != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.pressure))[index] = 0.0;
			}
		}
        }

	if(data.soundspeed != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.soundspeed))[index] = 0.0;
			}
		}
	}

	if(data.viscosity != NULL){
#pragma omp parallel for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.viscosity))[index] = 0.0;
			}
		}
	}

	if(data.xvel0 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.xvel0))[index] = 0.0;
			}
		}
	}

	if(data.xvel1 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.xvel1))[index] = 0.0;
			}
		}
	}

	if(data.yvel0 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.yvel0))[index] = 0.0;
			}
		}
	}

	if(data.yvel1 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.yvel1))[index] = 0.0;
			}
		}
	}



	if(data.work_array1 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array1))[index] = 0.0;
			}
		}
	}
	if(data.work_array2 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array2))[index] = 0.0;
			}
		}
        }
	if(data.work_array3 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array3))[index] = 0.0;
			}
		}
	}
	if(data.work_array4 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array4))[index] = 0.0;
			}
		}
	}
	if(data.work_array5 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array5))[index] = 0.0;
			}
		}
	}
	if(data.work_array6 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array6))[index] = 0.0;
			}
		}
	}
	if(data.work_array7 != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array7))[index] = 0.0;
			}
		}
	}


	if(data.vol_flux_x != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.vol_flux_x))[index] = 0.0;
			}
		}
	}
	if(data.vol_flux_y != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.vol_flux_y))[index] = 0.0;
			}
		}
	}

	if(data.mass_flux_x != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.mass_flux_x))[index] = 0.0;
			}
		}
	}
	if(data.mass_flux_y != NULL){
#pragma omp for private(k, j)
		for(k=y_min; k<=y_max; k++){
			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.mass_flux_y))[index] = 0.0;
			}
		}
	}
	
}

}

double sum_density0(struct data_obj data){
	double counter=0.0;
    int j, k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+2; k++){

		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
			counter+=(*(data.density0))[index];
		}
	}

	return counter;
}

double sum_energy0(struct data_obj data){
	double counter=0.0;
    int j,k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+2; k++){

		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
			counter+=(*(data.energy0))[index];
		}
	}

	return counter;
}

double sum_pressure(struct data_obj data){
	double counter=0.0;
    int j,k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+2; k++){

		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
			counter+=(*(data.pressure))[index];
		}
	}

	return counter;
}

double sum_soundspeed(struct data_obj data){
	double counter=0.0;
    int j,k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+2; k++){

		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
			counter+=(*(data.soundspeed))[index];
		}
	}

	return counter;
}

double sum_xvel1(struct data_obj data){
	double counter=0.0;
    int j,k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+3; k++){

		for(j=x_min-2; j<=x_max+3; j++){
			int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
			counter+=(*(data.xvel1))[index];
		}
	}

	return counter;
}

double sum_yvel1(struct data_obj data){
	double counter=0.0;
    int j,k;

    int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	for(k=y_min-2; k<=y_max+3; k++){

		for(j=x_min-2; j<=x_max+3; j++){
			int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
			counter+=(*(data.yvel1))[index];
		}
	}

	return counter;
}


double sum_2darray(struct data_obj data, double **arr, int x_bound, int y_bound){

	double counter=0.0;
	int j,k;

	int x_min, x_max, y_min, y_max;
	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);

	int x_max_b = x_max + x_bound;
	int y_max_b = y_max + y_bound;

	for(k=y_min-2; k<=y_max_b; k++){

		for(j=x_min-2; j<=x_max_b; j++){
			int index=FTNREF2D(j  ,k  ,x_max_b + 2,x_min-2,y_min-2);
			counter+=(*arr)[index];
		}
	}

	return counter;
}


void set_data(struct data_obj data){

    int x_min, x_max, y_min, y_max;
    int j,k;
    double dx,dy;

	printf("Setting Data\n");

	if((data.x_min == NULL) || (data.x_max == NULL) ||(data.y_min == NULL) ||(data.y_max == NULL) ){
		printf("Error: Setting with no problem size specified\n");
        return;
	}

	x_min=*(data.x_min);
	x_max=*(data.x_max);
	y_min=*(data.y_min);
	y_max=*(data.y_max);



	// Set the initial data

	dx=(10.0)/(x_max-x_min+1);
	dy=(10.0)/(y_max-y_min+1);

	if(data.vertexx != NULL){
		for(j=x_min-2; j<=x_max+3; j++){
			int index=FTNREF1D(j,x_min-2);
			(*(data.vertexx))[index] = 0.0+dx*((1.0*j)-x_min);
		}
	}

	if(data.vertexy != NULL){
		for(k=y_min-2; k<=y_max+3; k++){
			int index=FTNREF1D(k,y_min-2);
			(*(data.vertexy))[index] = 0.0+dy*((1.0*k)-y_min);
		}
	}


	if(data.cellx != NULL){
		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF1D(j,x_min-2);
			(*(data.cellx))[index] = 0.5*((*(data.vertexx))[index]+(*(data.vertexx))[index+1]);
		}

	}

	if(data.celly != NULL){
		for(k=y_min-2; k<=y_max+2; k++){
			int index=FTNREF1D(k,y_min-2);
			(*(data.celly))[index] = 0.5*((*(data.vertexy))[index]+(*(data.vertexy))[index+1]);
		}

	}


	if(data.vertexdx != NULL){
		for(j=x_min-2; j<=x_max+3; j++){
			int index=FTNREF1D(j,x_min-2);
			(*(data.vertexdx))[index] = dx;
		}

	}

	if(data.vertexdy != NULL){
		for(k=y_min-2; k<=y_max+3; k++){
			int index=FTNREF1D(k,y_min-2);
			(*(data.vertexdy))[index] = dy;
		}

	}


	if(data.celldx != NULL){
		for(j=x_min-2; j<=x_max+2; j++){
			int index=FTNREF1D(j,x_min-2);
			(*(data.celldx))[index] = dx;
		}

	}

	if(data.celldy != NULL){
		for(k=y_min-2; k<=y_max+2; k++){
			int index=FTNREF1D(k,y_min-2);
			(*(data.celldy))[index] = dy;
		}

	}


	if(data.xarea != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				int index2=FTNREF1D(k,y_min-2);

				(*(data.xarea))[index1] = (*(data.celldy))[index2];
			}
		}

	}

	if(data.yarea != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				int index2=FTNREF1D(j,x_min-2);

				(*(data.yarea))[index1] = (*(data.celldx))[index2];
			}
		}

	}

	if(data.volume != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);

				(*(data.volume))[index] = dx*dy;
			}
		}

	}



	if(data.density0 != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);

				double radius=sqrt( pow((j*1.0)*dx-5.0,2.0) + pow((k*1.0)*dy-5.0,2.0) );
				if(radius<=2.5){
					(*(data.density0))[index] = 2.0-(radius*2.0/10.0);
				}else{
					(*(data.density0))[index] = 1.0-((radius-5.0)*1.0/20.0);
				}
			}

		}

	}

	if(data.density1 != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);

				double radius=sqrt( pow((j*1.0)*dx-5.0,2.0) + pow((k*1.0)*dy-5.0,2.0) );
				if(radius<=2.5){
					(*(data.density1))[index] = 2.0-(radius*2.0/10.0);
				}else{
					(*(data.density1))[index] = 1.0-((radius-5.0)*1.0/20.0);
				}
			}

		}

	}

	if(data.energy0 != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);

				double radius=sqrt( pow((j*1.0)*dx-5.0,2.0) + pow((k*1.0)*dy-5.0,2.0) );
				if(radius<=2.5){
					(*(data.energy0))[index] = 2.0-(radius*2.0/10.0);
				}else{
					(*(data.energy0))[index] = 1.0-((radius-5.0)*1.0/20.0);
				}
			}

		}

	}

	if(data.energy1 != NULL){

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);

				double radius=sqrt( pow((j*1.0)*dx-5.0,2.0) + pow((k*1.0)*dy-5.0,2.0) );
				if(radius<=2.5){
					(*(data.energy1))[index] = 2.0-(radius*2.0/10.0);
				}else{
					(*(data.energy1))[index] = 1.0-((radius-5.0)*1.0/20.0);
				}
			}

		}

	}

	if(data.pressure != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				(*(data.pressure))[index] = (1.4-1.0)*(*(data.density0))[index]*(*(data.energy0))[index];
			}

		}
	}

	if(data.soundspeed != NULL){
		//*(data.soundspeed) = (double *) malloc(size_small_3d);
		/*
		 * v=1.0_8/density0(j,k)
        pressurebyenergy=(1.4_8-1.0_8)*density0(j,k)
        pressurebyvolume=-density0(j,k)*pressure(j,k)
        sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume)
        soundspeed(j,k)=SQRT(sound_speed_squared)
		 */

		for(k=y_min-2; k<=y_max+2; k++){

			for(j=x_min-2; j<=x_max+2; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				double v = 1.0/(*(data.density0))[index];
				double pressurebyenergy = (1.4-1.0)*(*(data.density0))[index];
				double pressurebyvolume=-(*(data.density0))[index]*(*(data.pressure))[index];
				double sound_speed_squared=v*v*((*(data.pressure))[index]*pressurebyenergy-pressurebyvolume);

				(*(data.soundspeed))[index] = sqrt(sound_speed_squared);
			}

		}
	}

	if(data.dt != NULL){

	    double dt = 0.0;
		(*(data.dt))=0.0;
		double width = MIN(dx,dy);
		for(k=y_min; k<=y_max; k++){

			for(j=x_min; j<=x_max; j++){
				int index=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				if((*(data.soundspeed))[index]>dt) dt = (*(data.soundspeed))[index];
			}
		}
		(*(data.dt))=width*0.7/dt;

	}

	if(data.xvel0 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);

				double x = (1.0*j)*dx-5.0;
				double y = (1.0*k)*dy-5.0;

				double radius=sqrt(pow(x,2)+pow(y,2));
				double mult = 0.0;
				if(x<=0.0 && y <= 0.0 ) mult = -1.0;
				if(x<=0.0 && y > 0.0 ) mult = -1.0;
				if(x>0.0 && y <= 0.0 ) mult = 1.0;
				if(x>0.0 && y > 0.0 ) mult = 1.0;

				double theta=0.0;
				if(x!=0.0){
					theta=atan(y/x);
				}else{
					theta=atan(y/-0.000000001);
				}

				if(radius<=2.5){
					(*(data.xvel0))[index] = mult*(2.0-(radius*2.0/10.0))*sin(theta);
				}else{
					(*(data.xvel0))[index] = mult*(1.0-((radius-5.0)*1.0/20))*sin(theta);
				}

			}

		}
	}

	if(data.yvel0 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);

				double x = (1.0*j)*dx-5.0;
				double y = (1.0*k)*dy-5.0;

				double radius=sqrt(pow(x,2)+pow(y,2));
				double mult = 0.0;
				if(x<=0.0 && y <= 0.0 ) mult = -1.0;
				if(x<=0.0 && y > 0.0 ) mult = -1.0;
				if(x>0.0 && y <= 0.0 ) mult = 1.0;
				if(x>0.0 && y > 0.0 ) mult = 1.0;

				double theta=0.0;
				if(x!=0.0){
					theta=atan(y/x);
				}else{
					theta=atan(y/-0.000000001);
				}
				radius = sqrt(pow((j*1.0)*dx-5.0,2.0)+pow((k*1.0)*dy-5.0,2.0));
				if(radius<=2.5){
					(*(data.yvel0))[index] = mult*(2.0-(radius*2.0/10.0))*cos(theta);
				}else{
					(*(data.yvel0))[index] = mult*(1.0-((radius-5.0)*1.0/20))*cos(theta);
				}

			}

		}
	}

	if(data.xvel1 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);

				double x = (1.0*j)*dx-5.0;
				double y = (1.0*k)*dy-5.0;

				double radius=sqrt(pow(x,2)+pow(y,2));
				double mult = 0.0;
				if(x<=0.0 && y <= 0.0 ) mult = -1.0;
				if(x<=0.0 && y > 0.0 ) mult = -1.0;
				if(x>0.0 && y <= 0.0 ) mult = 1.0;
				if(x>0.0 && y > 0.0 ) mult = 1.0;

				double theta=0.0;
				if(x!=0.0){
					theta=atan(y/x);
				}else{
					theta=atan(y/-0.000000001);
				}

				if(radius<=2.5){
					(*(data.xvel1))[index] = mult*(2.0-(radius*2.0/10.0))*sin(theta);
				}else{
					(*(data.xvel1))[index] = mult*(1.0-((radius-5.0)*1.0/20.0))*sin(theta);
				}

			}

		}
	}

	if(data.yvel1 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

					for(j=x_min-2; j<=x_max+3; j++){
						int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);

						double x = (1.0*j)*dx-5.0;
						double y = (1.0*k)*dy-5.0;

						double radius=sqrt(pow(x,2)+pow(y,2));
						double mult = 0.0;
						if(x<=0.0 && y <= 0.0 ) mult = -1.0;
						if(x<=0.0 && y > 0.0 ) mult = -1.0;
						if(x>0.0 && y <= 0.0 ) mult = 1.0;
						if(x>0.0 && y > 0.0 ) mult = 1.0;

						double theta=0.0;
						if(x!=0.0){
							theta=atan(y/x);
						}else{
							theta=atan(y/-0.000000001);
						}
						if(radius<=2.5){
							(*(data.yvel1))[index] = mult*(2.0-(radius*2.0/10.0))*cos(theta);
						}else{
							(*(data.yvel1))[index] = mult*(1.0-((radius-5.0)*1.0/20.0))*cos(theta);
						}

					}

				}
	}


	if(data.viscosity != NULL){


		for(k=y_min; k<=y_max; k++){

			for(j=x_min; j<=x_max; j++){
				int index1=FTNREF2D(j  ,k    ,x_max+4,x_min-2,y_min-2);
				int index1p2=FTNREF2D(j  ,k    ,x_max+5,x_min-2,y_min-2);
				int index2=FTNREF2D(j+1  ,k  ,x_max+4,x_min-2,y_min-2);
				int index2p2=FTNREF2D(j+1  ,k  ,x_max+5,x_min-2,y_min-2);
				int index3=FTNREF2D(j  ,k+1  ,x_max+4,x_min-2,y_min-2);
				int index3p2=FTNREF2D(j  ,k+1  ,x_max+5,x_min-2,y_min-2);
				int index4=FTNREF2D(j+1  ,k+1,x_max+4,x_min-2,y_min-2);
				int index4p2=FTNREF2D(j+1  ,k+1,x_max+5,x_min-2,y_min-2);
				int index5=FTNREF2D(j-1  ,k  ,x_max+4,x_min-2,y_min-2);
				int index5p2=FTNREF2D(j-1  ,k  ,x_max+5,x_min-2,y_min-2);
				int index6=FTNREF2D(j  ,k-1  ,x_max+4,x_min-2,y_min-2);
				int index6p2=FTNREF2D(j  ,k-1  ,x_max+5,x_min-2,y_min-2);


				int index7=FTNREF1D(j,x_min-2);

				int index8=FTNREF1D(k,y_min-2);

				double ugrad = ((*(data.xvel0))[index2p2]+(*(data.xvel0))[index4p2])-((*(data.xvel0))[index1p2]+(*(data.xvel0))[index3p2]);
				double vgrad = ((*(data.yvel0))[index3p2]+(*(data.yvel0))[index4p2])-((*(data.yvel0))[index1p2]+(*(data.yvel0))[index2p2]);

				double div = (*(data.celldx))[index7]*ugrad + (*(data.celldy))[index8]*vgrad;
				double strain2 = 0.5*((*(data.xvel0))[index3p2]+(*(data.xvel0))[index4p2]-(*(data.xvel0))[index1p2]-(*(data.xvel0))[index2p2])/ (*(data.celldy))[index8]\
						+ 0.5*((*(data.yvel0))[index2p2]+(*(data.yvel0))[index4p2]-(*(data.yvel0))[index1p2]-(*(data.yvel0))[index3p2])/ (*(data.celldx))[index7];

				double pgradx=((*(data.pressure))[index2] - (*(data.pressure))[index5])/((*(data.celldx))[index7] + (*(data.celldx))[index7+1]);
				double pgrady=((*(data.pressure))[index3] - (*(data.pressure))[index6])/((*(data.celldy))[index8] + (*(data.celldy))[index8+1]);


				double pgradx2 = pgradx*pgradx;
				double pgrady2 = pgrady*pgrady;

				double limiter = ((0.5*ugrad/(*(data.celldx))[index7])*pgradx2+(0.5*vgrad/(*(data.celldy))[index8])*pgrady2 + strain2*pgradx*pgrady) \
						/MAX(pgradx2+pgrady2,1.0e-16);

				if((limiter>0.0)||(div>=0.0)){
					(*(data.viscosity))[index1] = 0.0;
				}else{
					pgradx = SIGN(MAX(1.0e-16,fabs(pgradx)),pgradx);
					pgrady = SIGN(MAX(1.0e-16,fabs(pgrady)),pgrady);
					double pgrad = sqrt(pgradx*pgradx+pgrady*pgrady);
					double xgrad = fabs((*(data.celldx))[index7]*pgrad/pgradx);
					double ygrad = fabs((*(data.celldy))[index8]*pgrad/pgrady);
					double grad = MIN(xgrad,ygrad);
					double grad2 = grad*grad;

					(*(data.viscosity))[index1] = 2.0*(*(data.density0))[index1]*grad2*limiter*limiter;

				}


			}
		}

	}


	if(data.work_array1 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array1))[index] = 0.0;
			}

		}
	}
	if(data.work_array2 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array2))[index] = 0.0;
			}

		}
	}
	if(data.work_array3 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array3))[index] = 0.0;
			}

		}
	}
	if(data.work_array4 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array4))[index] = 0.0;
			}

		}
	}
	if(data.work_array5 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array5))[index] = 0.0;
			}

		}
	}
	if(data.work_array6 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array6))[index] = 0.0;
			}

		}
	}
	if(data.work_array7 != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min-2; k<=y_max+3; k++){

			for(j=x_min-2; j<=x_max+3; j++){
				int index=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				(*(data.work_array7))[index] = 0.0;
			}

		}
	}


	if(data.vol_flux_x != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min; k<=y_max; k++){

			for(j=x_min; j<=x_max+1; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				int index2=FTNREF2D(j  ,k+1  ,x_max+5,x_min-2,y_min-2);

				(*(data.vol_flux_x))[index1] = 0.25*(*(data.dt))*(*(data.xarea))[index1] *( (*(data.xvel0))[index1] + (*(data.xvel0))[index2] + (*(data.xvel1))[index1] + (*(data.xvel1))[index2]) ;
			}

		}
	}


	if(data.vol_flux_y != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min; k<=y_max+1; k++){

			for(j=x_min; j<=x_max; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				int index1p2=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				int index2=FTNREF2D(j+1  ,k  ,x_max+5,x_min-2,y_min-2);

				(*(data.vol_flux_y))[index1] = 0.25*(*(data.dt))*(*(data.yarea))[index1] *( (*(data.yvel0))[index1p2] + (*(data.yvel0))[index2] + (*(data.yvel1))[index1p2] + (*(data.yvel1))[index2]) ;
			}

		}
	}


	if(data.vol_flux_x != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min; k<=y_max; k++){

			for(j=x_min; j<=x_max+1; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2);
				int index2=FTNREF2D(j-1  ,k  ,x_max+4,x_min-2,y_min-2);

				(*(data.mass_flux_x))[index1] = (*(data.vol_flux_x))[index1] * (*(data.density1))[index2];
			}

		}
	}

	if(data.vol_flux_y != NULL){
		//*(data.pressure) = (double *) malloc(size_small_3d);

		for(k=y_min; k<=y_max+1; k++){

			for(j=x_min; j<=x_max; j++){
				int index1=FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2);
				int index2=FTNREF2D(j  ,k-1  ,x_max+4,x_min-2,y_min-2);

				(*(data.mass_flux_y))[index1] = (*(data.vol_flux_y))[index1] * (*(data.density1))[index2];
			}

		}
	}

}


