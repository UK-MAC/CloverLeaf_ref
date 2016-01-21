/*
 * data.h
 *
 *  Created on: 18 Jan 2016
 *      Author: ofjp
 */

#ifndef DATA_H_
#define DATA_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ftocmacros.h"

struct data_obj{
	/*
	  INTEGER :: x_min,x_max,y_min,y_max
	  REAL(KIND=8),OPTIONAL :: dt
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: cellx(:),celly(:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vertexdx(:),vertexdy(:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: celldx(:),celldy(:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xarea(:,:),yarea(:,:),volume(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: density0(:,:),density1(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: energy0(:,:),energy1(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: pressure(:,:),viscosity(:,:),soundspeed(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xvel0(:,:),yvel0(:,:),xvel1(:,:),yvel1(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_x(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_y(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_x(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_y(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array1(:,:),work_array2(:,:),work_array3(:,:),work_array4(:,:)
	  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array5(:,:),work_array6(:,:),work_array7(:,:)

	 */

	int *x_min;
	int *x_max;
	int *y_min;
	int *y_max;
	double *dt;
	double **cellx;
	double **celly;
	double **vertexdx;
	double **vertexdy;
	double **celldx;
	double **celldy;
	double **xarea;
	double **yarea;
	double **volume;
	double **density0;
	double **density1;
	double **energy0;
	double **energy1;
	double **pressure;
	double **viscosity;
	double **soundspeed;
	double **xvel0;
	double **yvel0;
	double **xvel1;
	double **yvel1;
	double **vol_flux_x;
	double **vol_flux_y;
	double **mass_flux_x;
	double **mass_flux_y;
	double **work_array1;
	double **work_array2;
	double **work_array3;
	double **work_array4;
	double **work_array5;
	double **work_array6;
	double **work_array7;

	double **vertexx;
	double **vertexy;


};


void init_data(struct data_obj *data);


void allocate_data(struct data_obj data);


double sum_density0(struct data_obj data);

double sum_energy0(struct data_obj data);

double sum_pressure(struct data_obj data);

double sum_soundspeed(struct data_obj data);

double sum_xvel1(struct data_obj data);

double sum_yvel1(struct data_obj data);


double sum_2darray(struct data_obj data, double **arr, int x_bound, int y_bound);


void set_data(struct data_obj data);



#endif /* DATA_H_ */
