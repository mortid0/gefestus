#include "initialStates.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "params.h"
#include "utils.h"
#include "pictures.h"

int initialize_walls();
  
//read data from data_file
void read_data()
{	
	FILE *fp;
	char str[8024], con[32], reg_exp[128];
	con[0]='\0';
	int i;
	if ((fp = fopen(data_file, "r"))==NULL)
	{
		printf("Can't open file %s\n",data_file);
		exit(1);
	}
	while (!feof(fp))
	{
		fscanf(fp, "%s",con);
		strcat(con," ");
		strcat(str, con);
	}
	fclose(fp);
	char res[N_MATCHES][MATCH_LEN];//number of matches
	match(res, "krp\\ *=\\ *([0-9]+)", str);
	KRP = atoi(res[1]);
	match(res, "kx\\ *=\\ *([0-9]+)", str);
	KX = atoi(res[1]);
	match(res, "ky\\ *=\\ *([0-9]+)", str);
	KY = atoi(res[1]);
	match(res, "nm\\ *=\\ *([0-9]+)", str);
	NM = atoi(res[1]);
	match(res, "nc\\ *=\\ *([0-9]+)", str);
	NC = atoi(res[1]);
	printf("krp=%i\nnc=%i\nnm=%i\nkx=%i\nky=%i\n",KRP,NC,NM,KX,KY);
	conc = malloc(NC*sizeof(double));
	temp = malloc(NC*sizeof(double));
	diam = malloc(NC*sizeof(double));
	mass = malloc(NC*sizeof(double));
	num_rot_degrees = malloc(NC*sizeof(int));
	num_vib_degrees = malloc(NC*sizeof(int));
 	ctmp = malloc(KRP*sizeof(double));
 	cdnn = malloc(KRP*sizeof(double));
	for (i = 0; i < NC; i++)
	{
		sprintf(reg_exp, "conc%d\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", (i+1));
		match(res, reg_exp, str);
		conc[i] = atof(res[1]);
		sprintf(reg_exp, "temp%d\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", (i+1));
		match(res, reg_exp, str);
		temp[i] = atof(res[1]);
		sprintf(reg_exp, "diam%d\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", (i+1));
		match(res, reg_exp, str);
		diam[i] = atof(res[1]);
		sprintf(reg_exp, "mass%d\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", (i+1));
		match(res, reg_exp, str);
		mass[i] = atof(res[1]);
		sprintf(reg_exp, "rot_deg%d\\ *=\\ *([0-9]+)", (i+1));
		match(res, reg_exp, str);
		num_rot_degrees[i] = atoi(res[1]);
		sprintf(reg_exp, "vib_deg%d\\ *=\\ *([0-9]+)", (i+1));
		match(res, reg_exp, str);
		num_vib_degrees[i] = atoi(res[1]);
	}
	match(res, "IS_INJECT\\ *=\\ *([0-9]+)", str);
	IS_INJECT = atoi(res[1]);
      //if (0 == IS_INJECT) return;
	match(res, "CONC_INJ\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", str);
	CONC_INJ = atof(res[1]);
	match(res, "TEMP_INJ\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", str);
	TEMP_INJ = atof(res[1]);
	match(res, "MASS_INJ\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", str);
	MASS_INJ = atof(res[1]);
	match(res, "DIAM_INJ\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", str);
	DIAM_INJ = atof(res[1]);
	match(res, "TIME_INJ\\ *=\\ *([0-9]+\\.*[0-9]*E?-?[0-9]+)", str);
	TIME_INJ = atof(res[1]);
	match(res, "CELL_INJ\\ *=\\ *([0-9]+)", str);
	CELL_INJ = atoi(res[1]);
	match(res, "dtp\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	dtp = atof(res[1]);
	match(res, "EXPERIMENT_TIME\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	EXPERIMENT_TIME = atof(res[1]);
	match(res, "SIZEX\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	SIZEX = atof(res[1]);
	match(res, "SIZEY\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	SIZEY = atof(res[1]);
	match(res, "SIZEZ\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	SIZEZ = atof(res[1]);
      /*	match(res, "kCollision\\ *=\\ *([0-9]+\\.*[0-9]*)", str);
	kCollision = atof(res[1]);*/
	match(res, "TEMP_FREQ\\ *=\\ *([0-9]+)", str);
	SHOW_TEMP_FREQ = atoi(res[1]);
	match(res, "CONC_FREQ\\ *=\\ *([0-9]+)", str);
	SHOW_CONC_FREQ = atoi(res[1]);
	match(res, "IMPULSE_FREQ\\ *=\\ *([0-9]+)", str);
	SHOW_IMPULSE_FREQ = atoi(res[1]);
	match(res, "MASS_FREQ\\ *=\\ *([0-9]+)", str);
	SHOW_MASS_FREQ = atoi(res[1]);
	  
//	printf("is_inj %i\nconc %f\ntemp %f\nmass %f\n diam %f\ntime %f\ncell %i\n", IS_INJECT, CONC_INJ, TEMP_INJ, MASS_INJ, DIAM_INJ, TIME_INJ, CELL_INJ);
	  
}
  
  // ��������� ������, ��� 5 ������ ������ ������
  // ��� �� ��������� buf -�������� ������, dtp_var, ������� �������� ������(left/right/up/down_wall_nalip), ������������� ��������� kNalip
  int create_camera()
{
	int krp, nc, kx, ky;
	camera = malloc(KRP*sizeof(PARTICLE ****));
	for (krp=0; krp<KRP; krp++)
	{
		camera[krp] = malloc(KX*sizeof(PARTICLE ***));
		for (kx=0; kx<KX; kx++)
		{
			camera[krp][kx] = malloc(KY*sizeof(PARTICLE **));
			for (ky = 0; ky<KY; ky++)
			{
				camera[krp][kx][ky] = malloc(NC*sizeof(PARTICLE *));
				for (nc=0; nc<NC; nc++)
				{
					camera[krp][kx][ky][nc] = malloc(NM*sizeof(PARTICLE));
				}
			}
		}
	}
	kNalip = malloc(NC*sizeof(double));
	svib = malloc(NC*sizeof(double));
	dtp_var = malloc(KRP*sizeof(double));
	temp_plate = malloc(KX*sizeof(double *));
	buf = malloc(KX*sizeof(PARTICLE **));
	buf_cell = malloc(NC*sizeof(PARTICLE *));
	for (nc=0; nc < NC; nc++)
	{
		buf_cell[nc] = malloc(NM*sizeof(PARTICLE));
	}
	force_field = malloc(KX * sizeof(VECTOR *));
	avr_vx = malloc(NC * sizeof(double));
	avr_vy = malloc(NC * sizeof(double));
	for (kx=0; kx<KX; kx++)
	{
		force_field[kx] = malloc( KY * sizeof(VECTOR));
		buf[kx] = malloc(KY*sizeof(PARTICLE *));
		for (ky=0; ky<KY; ky++)
		{
			buf[kx][ky] = malloc(NM*sizeof(PARTICLE));
		}
		temp_plate[kx] = malloc(KY*sizeof(double));
	}
	left_wall_nalip = malloc(KY*sizeof(double **));
	right_wall_nalip = malloc(KY*sizeof(double **));
	for (ky = 0; ky < KY; ky++)
	{
		left_wall_nalip[ky] = malloc(NC*sizeof(double *));
		right_wall_nalip[ky] = malloc(NC*sizeof(double *));
		for (nc = 0; nc < NC; nc++)
		{
			left_wall_nalip[ky][nc] = malloc(NALIP_MASS_COUNT*sizeof(double));
			right_wall_nalip[ky][nc] = malloc(NALIP_MASS_COUNT*sizeof(double));
		}
	}
	up_wall_nalip = malloc(KX*sizeof(double **));
	down_wall_nalip = malloc(KX*sizeof(double **));
	for (kx = 0; kx < KX; kx++)
	{
		up_wall_nalip[kx] = malloc(NC*sizeof(double *));
		down_wall_nalip[kx] = malloc(NC*sizeof(double *));
		for (nc = 0; nc < NC; nc++)
		{
			up_wall_nalip[kx][nc] = malloc(NALIP_MASS_COUNT*sizeof(double));
			down_wall_nalip[kx][nc] = malloc(NALIP_MASS_COUNT*sizeof(double));
		}
	}
	temp_field = malloc(NC * sizeof(double **));
	conc_field = malloc(NC * sizeof(double **));
	mass_field = malloc(NC * sizeof(double **));
	impulse_field = malloc(NC * sizeof(double **));
	distr_x = malloc(NC * sizeof(double *));
	distr_y = malloc(NC * sizeof(double *));
	for (nc = 0; nc < NC; nc++)
	{
		distr_x[nc] = malloc(MAX_PLOT_VEL * sizeof(double));
		distr_y[nc] = malloc(MAX_PLOT_VEL * sizeof(double));
		temp_field[nc] = malloc(KX * sizeof(double *));
		conc_field[nc] = malloc(KX * sizeof(double *));
		mass_field[nc] = malloc(KX * sizeof(double *));
		impulse_field[nc] = malloc(KX * sizeof(double *));
		for (kx = 0; kx < KX; kx++)
		{
				temp_field[nc][kx] = malloc(KY * sizeof(double));
				conc_field[nc][kx] = malloc(KY * sizeof(double));
				mass_field[nc][kx] = malloc(KY * sizeof(double));
				impulse_field[nc][kx] = malloc(KY * sizeof(double));
		}
	}
	return 0;
}
  
  //initialize camera with cam values and maxwell distribution
  int init()
{
	double d[NC];
	int i,j,k,l,m;
	for ( i = 0; i < NC; i++)
	{
		printf("conc#%i = %e\nmass#%i = %e\ndiam#%i = %e\nvmp#%i = %e\ntemp#%i = %e\nnum rot deg#%i = %i\nnum vib deg#%i = %i\n", i, conc[i], i, mass[i], i, diam[i], i, d[i], i, temp[i], i, num_rot_degrees[i], i, num_vib_degrees[i]);
		d[i] = sqrt(2.0*K*temp[i]/mass[i]); 
		
	}
	for (i = 0; i < 48; i++)
	{
		elev[0][i] = 3395*(i+0.5 - 0.006126*pow(i+0.5,2) + 0.00000318*pow(i+0.5,3)) - 1692.3;
		printf("elev[0][%i] = %e;\n",i,elev[0][i]);
	}
	for (k = 0; k<KRP; k++) 
	{
		for (i = 0; i < KX; i++) 
		{
			for (j = 0; j < KY; j++) 
			{
				for (l = 0; l<NC; l++) 
				{
					double conc_rem = conc[l], curr_conc;
					for (m = 0; m<NM; m++) 
					{
			  //camera[k][i][j][l][m].reflect_field = 0;
						camera[k][i][j][l][m].rot_degr = num_rot_degrees[l];
						camera[k][i][j][l][m].rot_energy = -log(frand()) * K * temp[l];
						camera[k][i][j][l][m].mass = mass[l];
						camera[k][i][j][l][m].diam = diam[l];
						if (NM-1 == m) {camera[k][i][j][l][m].weight = conc_rem;}
						else {
							curr_conc = 2.0*frand()*conc[l]/NM;
							if (curr_conc > conc_rem) {curr_conc = conc_rem / 2.0;}
							camera[k][i][j][l][m].weight = curr_conc;
							conc_rem -= curr_conc;
						}
						camera[k][i][j][l][m].v.x = get_hi2(0, d[l]);
						camera[k][i][j][l][m].v.y = sqrt(-log(frand()))*sin(2.0*PI*frand())*d[l];//get_hi2(0, d[l]);
						camera[k][i][j][l][m].v.z = sqrt(-log(frand()))*sin(2.0*PI*frand())*d[l];//get_hi2(0, d[l]);
			  
					}
				}
			}
		}
	}
      // 	printf("Den raio = %.5f\n", camera[0][0][0][0][0].weight/camera[0][0][0][0][1].weight);
	initialize_walls();
	avzv = 0;
	colzv = 0;
	return 0;
}

//initialize cam with values from file_name
int init_from_file(char *file_name)
{
	FILE *f;
	PARTICLE *p;
	int krp, nc, kx, ky, nm;
	f = fopen(file_name, "r");
	if (NULL == f)
	{
		printf("File %s not found!\n", file_name);
		return (-1);
	}
	while (!feof(f))
	{
		fscanf(f, "%i %i %i %i %i ", &krp, &kx, &ky, &nc, &nm);
		p = &camera[krp][kx][ky][nc][nm];
		fscanf(f, "%lf %lf %lf %lf %lf %lf\n", &(p->weight), &(p->mass), &(p->diam), &(p->v.x), &(p->v.y), &(p->v.z));
	}
	fclose(f);
	return 0;
}
 
//initialize camera with cam values and maxwell distribution, 1 particle move up, 0 - down
int init_2flows() 
{
	double d[NC];
	int i,j,k,l,m;
	for ( i = 0; i < NC; i++)
	{
		d[i] = sqrt(2.0*temp[i]/mass[i]); 
		printf("conc#%i = %e\nmass#%i = %e\ndiam#%i = %e\nvmp#%i = %e\ntemp#%i = %e\nnum rot deg#%i = %i\nnum vib deg#%i = %i\n", i, conc[i], i, mass[i], i, diam[i], i, d[i], i, temp[i], i, num_rot_degrees[i], i, num_vib_degrees[i]);
	}
	for (k = 0; k<KRP; k++) 
	{
		for (i = 0; i < KX; i++) 
		{
			for (j = 0; j < KY; j++) 
			{
				for (l = 0; l<NC; l++) 
				{
					for (m = 0; m<NM; m++) 
					{
			  //camera[k][i][j][l][m].reflect_field = 0;
						camera[k][i][j][l][m].rot_degr = num_rot_degrees[l];
						camera[k][i][j][l][m].rot_energy = 0;//-log(frand()) * K * temp[l];
						camera[k][i][j][l][m].mass = mass[l];
						camera[k][i][j][l][m].diam = diam[l];
						camera[k][i][j][l][m].weight = conc[l]/NM;
						
						if (0 == l)
						{
						    if ((i >= KX/2) && (j <= 4*KY/5) && (j>KY/3))
						    {
							camera[k][i][j][l][m].v.x = get_hi2(-temp[l],10);
							camera[k][i][j][l][m].v.y = get_hi2(0,10);
							camera[k][i][j][l][m].v.z = get_hi2(0,10);
							camera[k][i][j][l][m].weight = conc[l]/NM;
						    }
						    else {camera[k][i][j][l][m].weight = 0;/*
							camera[k][i][j][l][m].v.x = sin(2.0*PI*frand())*temp[l]/10;
							camera[k][i][j][l][m].v.y = sin(2.0*PI*frand())*temp[l]/10;
							camera[k][i][j][l][m].v.z = sin(2.0*PI*frand())*temp[l]/10;*/
							}
						}
						
						if (1 == l)
						{
						    if ((i < KX/2) && (j >= KY/5) && (j<2*KY/3))
						    {
							camera[k][i][j][l][m].v.x = get_hi2(temp[l],10);
							camera[k][i][j][l][m].v.y = get_hi2(0,10);
							camera[k][i][j][l][m].v.z = get_hi2(0,10);
							camera[k][i][j][l][m].weight = conc[l]/NM;
						    }
						    else {camera[k][i][j][l][m].weight = 0;
// 							camera[k][i][j][l][m].v.x = sin(2.0*PI*frand())*temp[l]/10;
// 							camera[k][i][j][l][m].v.y = sin(2.0*PI*frand())*temp[l]/10;
// 							camera[k][i][j][l][m].v.z = sin(2.0*PI*frand())*temp[l]/10;
							}
						}
						if (2 == l) {camera[k][i][j][l][m].weight = 0;}
			  
					}
				}
			}
		}
	}
	initialize_walls();
	avzv = 0;
	colzv = 0;
	return 0;
}
  
  //  ������������� ���� 
  int initialize_walls()
{
	left_wall.temp = 1;
	left_wall.kAcc = 1;
	left_wall.kFriction = 1;
	right_wall.temp = 1;
	right_wall.kAcc = 1;
	right_wall.kFriction = 1;
	up_wall.temp = 1;
	up_wall.kAcc = 1;
	up_wall.kFriction = 1;
	down_wall.temp = 1;
	down_wall.kAcc = 1;
	down_wall.kFriction = 1;
	return 0;
}
  
