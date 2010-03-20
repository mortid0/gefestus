#include <stdio.h>
#include "params.h"
#include "algoritm.h"
#include "move.h"
#include "pictures.h"
#include "utils.h"
#include "initialStates.h"
#include "collision.h"
#include <time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include "filesystem.h"
	

//TODO vstavit ispolzovanie conc_leave 

void prepare_root();

void clear_up();

void run_single();

int experiment_step(int step);


void prepare_root() 
{
	make_dirs();
}

void clear_up() 
{
	printf("Clear start\n");
	fclose(plot_file);
	int krp, nc, kx, ky;
	printf("Plot file closed\n");
	for (krp=0; krp<KRP; krp++)
	{
		for (kx=0; kx<KX; kx++)
		{
			for (ky = 0; ky<KY; ky++)
			{
				for (nc=0; nc<NC; nc++)
				{
					free(camera[krp][kx][ky][nc]);
				}
				free(camera[krp][kx][ky]);
			}
			free(camera[krp][kx]);
		}
		free(camera[krp]);
	}
	free(camera);
	printf("Camera cleaned\n");
	free(kNalip);
	free(dtp_var);
	free(conc);
	free(temp);
	free(mass);
	free(diam);
	printf("Macroprams cleaned\n");
	for (kx=0; kx<KX; kx++)
	{
		for (ky=0; ky<KY; ky++)
		{
			free(buf[kx][ky]);
		}
		free(buf[kx]);
		free(temp_plate[kx]);
	}
	free(temp_plate);
	free(buf);
	printf("Plates cleaned\n");
	for (ky = 0; ky < KY; ky++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			free(left_wall_nalip[ky][nc]);
			free(right_wall_nalip[ky][nc]);
		}
		free(left_wall_nalip[ky]);
		free(right_wall_nalip[ky]);
	}
	free(left_wall_nalip);
	free(right_wall_nalip);
	for (kx = 0; kx < KX; kx++)
	{
		
		for (nc = 0; nc < NC; nc++)
		{
			free(up_wall_nalip[kx][nc]);
			free(down_wall_nalip[kx][nc]);
		}
		free(up_wall_nalip[kx]);
		free(down_wall_nalip[kx]);
	}
	free(up_wall_nalip);
	free(down_wall_nalip);
// 	free(temp_field);
// 	free(conc_field);
// 	free(mass_field);
		
	printf("Nalip walls cleaned\n");
}

#define TOTAL_EXP 1
#define diam_0 6000
#define diam_step 500

void get_cam_mass_conc()
{
	int nc, krp, kx, ky, nm;
	for (nc = 0; nc < NC; nc++)
	{
		fin_m[nc] = 0;
		fin_c[nc] = 0;
		for (krp = 0; krp < KRP; krp++)
		{
			for (kx = 0; kx < KX; kx++)
			{
				for (ky = 0; ky < KY; ky++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						fin_c[nc] += camera[krp][kx][ky][nc][nm].weight;
						fin_m[nc] += camera[krp][kx][ky][nc][nm].weight*camera[krp][kx][ky][nc][nm].mass;
					}
				}
			}
		}
	}
}

void run_single()
{
	int i, krp, exp_num;
	time_t start, finish;
	start = time((time_t *)NULL);
	printf("start = %s\n", ctime(&start));
	read_data();
	prepare_root();
	create_camera();
// 	normirovka();
	
	fclose(plot_file);
	plot_file = fopen("temp2/plot0.dat", "w+");
	fprintf(plot_file, "# DIAM  \t FIN_C[0]  \t FIN_M[0]  \t FIN_C[1]  \t FIN_M[1]  \t FIN_C[2]  \t FIN_M[2]  \t TOT_C  \t  TOT_M\n");
	for (exp_num = 0; exp_num < TOTAL_EXP; exp_num++)
	{
		double curr_diam = diam_0 + exp_num * diam_step;
		diam[0] = curr_diam;
		diam[1] = curr_diam;
		diam[2] = curr_diam;
//		init_2flows();
		init_from_file("temp2/data.in");
		double new_exp_time = EXPERIMENT_TIME;//KX * SIZEX / temp[0];
		int time_steps = (int)(new_exp_time / dtp) ;//(int)(EXPERIMENT_TIME/dtp);
		printf("Diam = %4.4e; Dens singe aster = %4.4e kg m^-3\n", curr_diam, mass[0]/(0.5 * 0.333 * PI *curr_diam*curr_diam*curr_diam));
		printf("Time steps = %i;\n", time_steps);
		for (krp = 0; krp < KRP; krp++)
		{
			dtp_var[krp] = dtp;
		}
		for (i=0; i <= time_steps; i++)
		{
			experiment_step(i);
			printf("*************** TIME = %e ***** STEP = %i************\n", i*dtp, i);
		}
		get_cam_mass_conc();
		fprintf(plot_file, "%4.4e;\t", curr_diam);
		for (i = 0; i < NC; i++)
		{
			fprintf(plot_file, "%4.4e;\t%4.4e;\t", fin_c[i], fin_m[i]);
			printf("fin_m[%i] = %5.5e;\t fin_c[%i] = %5.5e;\n", i, fin_m[i], i, fin_c[i]);
		}
		fprintf(plot_file,"\t%4.4e;\t%4.4e;\n", fin_c[0]+fin_c[1]+fin_c[2], fin_m[0]+fin_m[1]+fin_m[2]);
		printf("Exp #%i over; %e; %.5f;\n", exp_num, cell_temp, alfa);
	}
	clear_up();
	finish = time((time_t *)NULL);
	start = finish - start;
	printf("finish = %s\n", ctime(&finish));
	printf("total = %i\n", start);
}


int experiment_step(int step)
{
	int krp;
	int nc;

//   	mixture();
	print_pictures(step);
	
	for (krp = 0; krp < KRP; krp++)
	{
//        	inelastic_collision(krp);
	    	force(krp);
		for (nc = 0; nc < NC; nc++)
		{
		    move_camera(nc, krp);
		}
	}
// 	inject_left(step*dtp);
// 	inject_right(step*dtp);
	return 0;
}


int main(int argc, char *argv[]) 
{
	run_single();
	return 0;
}
