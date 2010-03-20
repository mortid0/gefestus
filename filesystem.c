#include "filesystem.h"
#include "params.h"
#include "pictures.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>


void make_dir_with_components(char *str)
{
	mode_t mode;
	mode = S_IRWXU;
	mkdir(str, mode);
	char str_comp_dir[40];
	int i;
	for (i = 0; i < NC; i++)
	{
		sprintf(str_comp_dir, "%s/%i", str, i);
		mkdir(str_comp_dir, mode);
	}
}

void make_dirs()
{
	if (SHOW_TEMP_FREQ > 0) make_dir_with_components(TEMP_DIR);
	if (SHOW_CONC_FREQ > 0) make_dir_with_components(CONC_DIR);
	make_dir_with_components("temp2/full_mass");
	make_dir_with_components("temp2/plotx");
	make_dir_with_components("temp2/ploty");
	if (SHOW_IMPULSE_FREQ > 0) make_dir_with_components(IMPULSE_DIR);
// 	if (SHOW_MASS_FREQ > 0) {int old_nc = NC; NC = 50; make_dir_with_components(MASS_DIR); NC = old_nc;}
	if (!(plot_file = fopen(PLOT_FILE, "w+")))
	{
		printf("#%i can't open file for write: %s\n", my_rank, PLOT_FILE);
		exit(1);
	}
}

void save_camera(int step)
{
	FILE *fp;
	char file_path[40];
	PARTICLE *p;
	int krp, nc, nm, x, y, nmc;
	sprintf(file_path, "%s/%i/%i.txt", SAVE_DIR, my_rank, step);
	if (!(fp = fopen(file_path, "w+")))
	{
		printf("#%i can't open file for write: %s\n", my_rank, file_path);
		exit(1);
	}
	for (krp = 0; krp < KRP; krp++)
	{
		for (x = 0; x < KX; x++)
		{
			for (y = 0; y < KY; y++)
			{
				for (nc = 0; nc < NC; nc++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						p = &camera[krp][x][y][nc][nm];
						fprintf(fp, "%f %f %f %f %f %f\n", p->mass, p->diam, p->weight, p->v.x, p->v.y, p->v.z);
					}
				}
			}
		}
	}
/*  for (krp = 0; krp < KRP; krp++)
	{
	fprintf(fp, "%f\n", dtp_var[krp]);
}*/
	for (nc = 0; nc < NC; nc++)
	{
		fprintf(fp, "%f\n", kNalip[nc]);
		for  (nmc = 0; nmc < NALIP_MASS_COUNT; nmc++)
		{
			for (x = 0; x < KX; x++)
			{
				fprintf(fp, "%f\n", up_wall_nalip[x][nc][nmc]);
				fprintf(fp, "%f\n", down_wall_nalip[x][nc][nmc]);
			}
			for (y = 0; y < KY; y++)
			{
				fprintf(fp, "%f\n", left_wall_nalip[y][nc][nmc]);
				fprintf(fp, "%f\n", right_wall_nalip[y][nc][nmc]);
			}
		}
	}
	fprintf(fp, "%f\n", conc_leave);
	fclose(fp);
}

void load_camera(int step)
{

	FILE *fp;
	char file_path[40];
	PARTICLE *p;
	int krp, nc, nm, x, y, nmc;
	sprintf(file_path, "%s/%i/%i.txt", SAVE_DIR, my_rank, step);
	if (!(fp = fopen(file_path, "r")))
	{
		printf("#%i can't open file for read: %s\n", my_rank, file_path);
		exit(1);
	}
	for (krp = 0; krp < KRP; krp++)
	{
		for (x = 0; x < KX; x++)
		{
			for (y = 0; y < KY; y++)
			{
				for (nc = 0; nc < NC; nc++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						p = &camera[krp][x][y][nc][nm];
						fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &p->mass, &p->diam, &p->weight, &p->v.x, &p->v.y, &p->v.z);
					}
				}
			}
		}
	}
/*  for (krp = 0; krp < KRP; krp++)
	{
	fscanf(fp, "%f\n", dtp_var[krp]);
}*/
	for (nc = 0; nc < NC; nc++)
	{
		fscanf(fp, "%lf\n", &kNalip[nc]);
		for  (nmc = 0; nmc < NALIP_MASS_COUNT; nmc++)
		{
			for (x = 0; x < KX; x++)
			{
				fscanf(fp, "%lf\n", &up_wall_nalip[x][nc][nmc]);
				fscanf(fp, "%lf\n", &down_wall_nalip[x][nc][nmc]);
			}
			for (y = 0; y < KY; y++)
			{
				fscanf(fp, "%lf\n", &left_wall_nalip[y][nc][nmc]);
				fscanf(fp, "%lf\n", &right_wall_nalip[y][nc][nmc]);
			}
		}
	}
	fscanf(fp, "%lf\n", &conc_leave);
	fclose(fp);
}
