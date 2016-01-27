/*
 * PdV_driver_c.c
 *
 *  Created on: 27 Jan 2016
 *      Author: ofjp
 */


#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[] ){

	int x_min, x_max, y_min, y_max;
	double dt;
	double *celldx;
	double *celldy;
	double *vertexdx;
	double *vertexdy;
	double *xarea;
	double *yarea;
	double *volume;
	double *density1;
	double *energy1;
	double *density0;
	double *pressure;
	double *soundspeed;
	double *energy0;
	double *xvel0;
	double *yvel0;
	double *xvel1;
	double *yvel1;
	double *vol_flux_x;
	double *vol_flux_y;
	double *mass_flux_x;
	double *mass_flux_y;
	double *post_ener;
	double *post_mass;
	double *pre_mass;
	double *advec_vol;
	double *ener_flux;
	double *pre_vol;
	double *post_vol;
	double *density_orig;
	double *energy_orig;

	double time_s;
	double time_e;
	double time_c;

	int iterations;
	int i;
	int sweep_number = 0;
	int direction = 0;
	int reset=0;

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

		if(strcmp(argv[i],"-reset") == 0){
			if(strcmp(argv[i+1],"on") == 0) reset=1;
			if(strcmp(argv[i+1],"off") == 0) reset=0;
		}
	}


	// Update Values

	x_min=1;
	y_min=1;
	x_max=x_size;
	y_max=y_size;

	// Print


	printf("PdV Kernel\n");
	printf("Mesh size %d %d\n", x_size, y_size);
	printf("Iterations %d\n", its);

	// Init Data Object


	struct data_obj Data;
	init_data(&Data);

	Data.x_min = &x_min;
	Data.x_max = &x_max;
	Data.y_min = &y_min;
	Data.y_max = &y_max;
	Data.dt = &dt;


	Data.celldx=&celldx;
	Data.celldy=&celldy;
	Data.vertexdx=&vertexdx;
	Data.vertexdy=&vertexdy;
	Data.xarea=&xarea;
	Data.yarea=&yarea;
	Data.volume=&volume;
	Data.density1=&density1;
	Data.energy1=&energy1;
	Data.density0=&density0;
	Data.pressure=&pressure;
	Data.soundspeed=&soundspeed;
	Data.energy0=&energy0;
	Data.xvel0=&xvel0;
	Data.yvel0=&yvel0;
	Data.xvel1=&xvel1;
	Data.yvel1=&yvel1;
	Data.vol_flux_x=&vol_flux_x;
	Data.vol_flux_y=&vol_flux_y;
	Data.mass_flux_x=&mass_flux_x;
	Data.mass_flux_y=&mass_flux_y;
	Data.work_array1=&post_ener;
	Data.work_array2=&post_mass;
	Data.work_array3=&pre_mass;
	Data.work_array4=&advec_vol;
	Data.work_array5=&ener_flux;
	Data.work_array6=&pre_vol;
	Data.work_array7=&post_vol;
	Data.reset_density=&density_orig;
	Data.reset_energy=&energy_orig;


	allocate_data(Data);
	set_data(Data);


	printf("Pre: Density: %f\n", sum_2darray(Data, Data.density1, 2, 2));
	printf("Pre: Energy: %f\n", sum_2darray(Data, Data.energy1, 2, 2));

	if(reset==1){
	    reset_store(Data, Data.density1, Data.reset_density);
	    reset_store(Data, Data.energy1, Data.reset_energy);
	}


	printf("Running Kernel\n");

	// Run Kernel
	time_c = 0.0;

	for(iterations=1;iterations<=its;iterations++){
		if(reset==1){
			reset_store(Data, Data.reset_density, Data.density1);
			reset_store(Data, Data.reset_energy, Data.energy1);
		}

		timer_c_(&time_s);

	    direction=1;
	    sweep_number=1;

	    advec_cell_kernel_c_(Data.x_min,
				Data.x_max,
				Data.y_min,
				Data.y_max,
				&direction,
				&sweep_number,
				*Data.vertexdx,
				*Data.vertexdy,
				*Data.volume,
				*Data.density1,
				*Data.energy1,
				*Data.mass_flux_x,
				*Data.vol_flux_x,
				*Data.mass_flux_y,
				*Data.vol_flux_y,
				*Data.work_array6,
				*Data.work_array7,
				*Data.work_array3,
				*Data.work_array2,
				*Data.work_array4,
				*Data.work_array1,
				*Data.work_array5);

	    direction=2;
	    sweep_number=2;

	    advec_cell_kernel_c_(Data.x_min,
				Data.x_max,
				Data.y_min,
				Data.y_max,
				&direction,
				&sweep_number,
				*Data.vertexdx,
				*Data.vertexdy,
				*Data.volume,
				*Data.density1,
				*Data.energy1,
				*Data.mass_flux_x,
				*Data.vol_flux_x,
				*Data.mass_flux_y,
				*Data.vol_flux_y,
				*Data.work_array6,
				*Data.work_array7,
				*Data.work_array3,
				*Data.work_array2,
				*Data.work_array4,
				*Data.work_array1,
				*Data.work_array5);
	timer_c_(&time_e);
	time_c += time_e-time_s;

	}

	// Print Result

	printf("Advec Cell time %f\n", time_c);
	printf("Density: %f\n", sum_2darray(Data, Data.density1, 2, 2));
	printf("Energy: %f\n", sum_2darray(Data, Data.energy1, 2, 2));


	// Clean up data
	free_data(Data);


	return 0;
}
