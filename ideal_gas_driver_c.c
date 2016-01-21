/*
 * ideal_gas_driver_c.c
 *
 *  Created on: 18 Jan 2016
 *      Author: ofjp
 */


#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[] ){

	int x_min, x_max, y_min, y_max;
	double *density0;
	double *energy0;
	double *pressure;
	double *soundspeed;

	double time_s;
	double time_e;

	int iterations;
	int i;

	// Defaults

	int x_size=100;
	int y_size=100;
	int its=1;


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
	}


	// Update Values

	x_min=1;
	y_min=1;
	x_max=x_size;
	y_max=y_size;

	// Print


	printf("Ideal Gas Kernel\n");
	printf("Mesh size %d %d\n", x_size, y_size);
	printf("Iterations %d\n", its);

	// Init Data Object


	struct data_obj Data;
	init_data(&Data);

	Data.x_min = &x_min;
	Data.x_max = &x_max;
	Data.y_min = &y_min;
	Data.y_max = &y_max;


	Data.density0 = &density0;
	Data.energy0 = &energy0;
	Data.pressure = &pressure;
	Data.soundspeed = &soundspeed;

	allocate_data(Data);
	set_data(Data);
	printf("Pre: Density: %f\n", sum_density0(Data));
	printf("Pre: Energy: %f\n", sum_energy0(Data));
	printf("Pre: Pressure: %f\n", sum_pressure(Data));
	printf("Pre: Soundspeed: %f\n", sum_soundspeed(Data));


	printf("Allocated Data - End\n");

	// Run Kernel
	timer_c_(&time_s);
	for(iterations=1;iterations<=its;iterations++){

		ideal_gas_kernel_c_(Data.x_min,Data.x_max,Data.y_min, Data.y_max, *Data.density0, *Data.energy0, *Data.pressure, *Data.soundspeed );

	}
	timer_c_(&time_e);

	// Print Result

	printf("Ideal gas time %f\n", time_e-time_s);
	printf("Density: %f\n", sum_density0(Data));
	printf("Energy: %f\n", sum_energy0(Data));
	printf("Pressure: %f\n", sum_pressure(Data));
	printf("Soundspeed: %f\n", sum_soundspeed(Data));


	// Clean up data
	free(density0);
	free(energy0);
	free(pressure);
	free(soundspeed);


	return 0;
}
