/*
 * accelerate_driver_c.c
 *
 *  Created on: 18 Jan 2016
 *      Author: ofjp
 */



#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[] ){

	int x_min, x_max, y_min, y_max;
	double *xarea;
	double *yarea;
	double *volume;
	double *celldx;
	double *celldy;
	double *density0;
	double *energy0;
	double *pressure;
	double *soundspeed;
	double *viscosity;
	double *xvel0;
	double *yvel0;
	double *xvel1;
	double *yvel1;
	double *work_array1;
	double *xvel_orig;
	double *yvel_orig;
	double *vertexx;
	double *vertexy;
	double dt;


	double time_s;
	double time_e;

	int iterations;
	int i;

	// Defaults

	int x_size=100;
	int y_size=100;
	int its=1;
	int reset=0;



    // Read Data
	for(i=1; i<argc; i++){
		if(strcmp(argv[i],"-nx") == 0){
			x_size=atoi(argv[i+1]);
		}

		if(strcmp(argv[i],"-ny") == 0){
			y_size=atoi(argv[i+1]);
		}

		if(strcmp(argv[i],"-its") == 0){
			its=atoi(argv[i+1]);
		}

		if(strcmp(argv[i],"-reset") == 0){
			if(strcmp(argv[i+1],"on") == 0){
				reset=1;
			}else{
				reset=0;
			}
		}
	}


	// Update Values

	x_min=1;
	y_min=1;
	x_max=x_size;
	y_max=y_size;

	// Print


	printf("Accelerate Kernel\n");
	printf("Mesh size %d %d\n", x_size, y_size);
	printf("Iterations %d\n", its);

	// Init Data Object


	struct data_obj Data;

	init_data(&Data);


	Data.dt = &dt;

	Data.x_min = &x_min;
	Data.x_max = &x_max;
	Data.y_min = &y_min;
	Data.y_max = &y_max;

	Data.xarea = &xarea;
	Data.yarea= &yarea;
	Data.volume= &volume;
	Data.celldx = &celldx;
	Data.celldy = &celldy;

	Data.density0 = &density0;
	Data.energy0 = &energy0;
	Data.pressure = &pressure;
	Data.soundspeed = &soundspeed;
	Data.viscosity = &viscosity;

	Data.xvel0 = &xvel0;
	Data.yvel0 = &yvel0;
	Data.xvel1 = &xvel1;
	Data.yvel1 = &yvel1;
	Data.work_array1 = &work_array1;

	Data.vertexx = &vertexx;
	Data.vertexy = &vertexy;


	allocate_data(Data);
	set_data(Data);


	printf("Running Kernel\n");

	// Run Kernel
	timer_c_(&time_s);
	for(iterations=1;iterations<=its;iterations++){

		//ideal_gas_kernel_c_(Data.x_min,Data.x_max,Data.y_min, Data.y_max, *Data.density0, *Data.energy0, *Data.pressure, *Data.soundspeed );
		accelerate_kernel_c_(Data.x_min,                  \
				Data.x_max,                  \
				Data. y_min,                  \
				Data.y_max,                  \
				Data.dt,                     \
				*Data.xarea,                  \
				*Data.yarea,                  \
				*Data.volume,                 \
				*Data.density0,               \
				*Data.pressure,               \
				*Data.viscosity,              \
				*Data.xvel0,                  \
				*Data.yvel0,                  \
				*Data.xvel1,                  \
				*Data.yvel1,                  \
				*Data.work_array1             );


	}
	timer_c_(&time_e);

	// Print Result

	printf("Accelerate time %f\n", time_e-time_s);
	printf("XVel1: %f\n", sum_2darray(Data, Data.xvel1, 3, 3));
	printf("YVel1: %f\n", sum_2darray(Data, Data.yvel1, 3, 3));


	// Clean up data
	free_data(Data);


	return 0;
}
