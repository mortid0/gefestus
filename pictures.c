#include "pictures.h"
#include "params.h"
#include <stdio.h>
#include <math.h>
#include "utils.h"



void show_field(int num, double (*cell_func)(int, int, int), char *func_dir,  int startx, int endx, int starty, int endy)
{
	FILE *fp;
	char file_name[60];
	int nc, x, y;
	double cell_res, sum = 0;
	for (nc = 0; nc < NC; nc++)
	{
		sum = 0;
		sprintf(file_name, "%s/%i/%i%s", func_dir, nc, num, ".txt");
		fp = fopen(file_name, "w+");
		for (y = starty; y <= endy; y++)
		{
			for (x = startx; x <= endx; x++)
			{
				cell_res = (*cell_func)(nc, x, y);
				sum += cell_res;
				fprintf(fp, "%.8e; %.8e; %.8e;\n", SIZEX*x, SIZEY * y,(cell_res));
			}
			fprintf(fp, "\n");
		}
		printf("#%i %s%d = %e, avr = %e;\n", my_rank, func_dir, nc, sum, sum / ((endx - startx + 1) * (endy - starty + 1)));
		/* 		fprintf(plot_file, "%.5f ", sum/(KX*KY)); */
		fclose(fp);
	}
	/* 	fprintf(plot_file, " %.5f\n", num*dtp); */
}

void show_temp_2flows(int num)
{
	FILE *fp;
	char file_name[60];
	int nc, x, y, krp, nm;
	double  norm;
	double tleft, tright, norml, normr;
	PARTICLE *p;
	VECTOR v0;
	nc = 0;

	sprintf(file_name, "%s/%i/%i%s", "temp2/temp2flows", nc, num, ".txt");
	fp = fopen(file_name, "w+");
	for (y = 0; y < KY ; y++)
	{
		for (x = 0; x < KX; x++)
		{
			norm = 0.0;
			v0.x = v0.y = v0.z = 0.0;
			norml = tleft = 0.0;
			normr = tright = 0.0;
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						v0.x += p->v.x * p->weight;
						v0.y += p->v.y * p->weight;
						v0.z += p->v.z * p->weight;
						norm += p->weight;
					}
				}
				multiply(&v0, &v0, (1.0/norm));
			}
/* 			sum_t = 0.0; */
			for (nc = 0; nc < NC; nc++)
			{
/* 				tleft = tright = 0; */
/* 				norml = normr = 0; */
				for (nm = 0; nm < NM; nm++)
				{
					/* 					printf("v0.x = %.5f\n", v0.x); */
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						if (p->v.x >= 0)
						{
							tleft += get_temp(p, &v0);
							norml += p->weight;
						}
						if (p->v.x < 0)
						{
							tright += get_temp(p, &v0);
							normr += p->weight;
						}
					}
				}
/* 				sum_t += (tleft+tright)/(norml+normr); */
			}
			fprintf(fp, "%.5f;   %.5f;          %.5f\n", (tleft / norml), (tright / normr), (tleft+tright)/(norml+normr));
		}
		fclose(fp);
	}
}

double get_freq(double temper)
{
    double rm = 0.5*mass[0];
    double cr = 2.0 * sqrt(2.0*K*temper/(rm*PI));
    double sigma = PI*diam[0]*diam[0];
    double fr = conc[0]*sigma*pow(cr, 2.0*viscos - 1)*(2.0/sqrt(PI))*gam(2.5 - viscos)*pow(2.0*K*temper/rm,1.0-viscos);
    return fr;
}

void show_temp_internal(int num)
{
	int nc, x, y, krp, nm;
	PARTICLE *p;
	double cs0[NC], csvs[NC][2], csv[NC], csr[NC];
	VECTOR cs[NC], cs2[NC];
	nc = 0;
	for (y = 0; y < KY ; y++)
	{
		for (x = 0; x < KX; x++)
		{
			for (krp = 0; krp < KRP; krp++)
			{
				cdnn[krp] = 0;
				ctmp[krp] = 0;
			}
			for (nc = 0; nc < NC; nc++)
			{
				cs0[nc] = 0;
				csvs[nc][1] = csvs[nc][0] = 0;
				csv[nc] = 0;csr[nc] = 0;
				cs[nc].x = cs[nc].y = cs[nc].z = 0;
				cs2[nc].x = cs2[nc].y = cs2[nc].z = 0;
				csr[nc] = 0;
				csv[nc] = 0;
				for (nm = 0; nm < NM; nm++)
				{
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						cs0[nc] += p->weight;
						cdnn[krp] += p->weight;
						ctmp[krp] += p->weight*p->mass*mag2(&(p->v));
						cs[nc].x += p->v.x * p->weight;
						cs[nc].y += p->v.y * p->weight;
						cs[nc].z += p->v.z * p->weight;
						cs2[nc].x += pow(p->v.x, 2) * p->weight;
						cs2[nc].y += pow(p->v.y, 2) * p->weight;
						cs2[nc].z += pow(p->v.z, 2) * p->weight;
						if (p->rot_degr > 0) {csr[nc] += p->rot_energy * p->weight;}
					}
				}
			}
			for (krp = 0; krp < KRP; krp++)
			{
				ctmp[krp] /= cdnn[krp]*3.0*K;
			}
			for (nc = 0; nc < NC; nc++)
			{
				cs0[nc] /= KRP;
				cs[nc].x  /= KRP;
				cs[nc].y /= KRP;
				cs[nc].z /= KRP;
				cs2[nc].x  /= KRP;
				cs2[nc].y /= KRP;
				cs2[nc].z /= KRP;
				csr[nc] /= KRP;
				csv[nc] /= KRP;
				csvs[nc][0]  /= KRP;
				csvs[nc][1] /= KRP;
			}
			alfa = cs0[1]/(2.0*cs0[0]+cs0[1]);
			double sn = 0, sm = 0, smcc = 0, sre = 0, srdf = 0;
			double tvib[NC];
			VECTOR smu;
			smu.x = smu.y = smu.z = 0.0;
			for (nc = 0; nc < NC; nc++)
			{
				sn += cs0[nc];
				sm += mass[nc]*cs0[nc];
				smu.x += mass[nc]*cs[nc].x; smu.y += mass[nc]*cs[nc].y; smu.z += mass[nc]*cs[nc].z;
				smcc += mass[nc]*(cs2[nc].x+cs2[nc].y+cs2[nc].z);
				sre += csr[nc];
				srdf += cs0[nc]*num_rot_degrees[nc];
				if (num_vib_degrees[nc] >0)
				{
					if (csvs[nc][1] > 0) 
					{
						tvib[nc] = REF_VIB_TEMP/log(csvs[nc][0]/csvs[nc][1]); 
						svib[nc] = 2.0 * csv[nc] / (cs0[nc] * K * tvib[nc]);
					} else
					{
						tvib[nc] = 1.0E-6;
						svib[nc] = 0;
					}
				}
			}
			double denn = sn;
			double den = denn * sm / sn;
			VECTOR vel, svel;
			vel.x = smu.x/sm;vel.y = smu.y/sm;vel.z = smu.z/sm;
			svel.x = vel.x; svel.y = vel.y; svel.z = vel.z;
			double uu = pow(vel.x,2) + pow(vel.y,2) + pow(vel.z,2);
			double tt = (smcc-sm*uu)/(3.0*K*sn);
			double trot = (2.0 / K) * sre / srdf;
			double svt = 0;
			double svdf = 1.0E-8;
			for (nc = 0; nc < NC; nc++)
			{
				if (num_vib_degrees[nc] > 0)
				{
					svt += tvib[nc]*svib[nc]*cs0[nc];
					svdf += svib[nc]*cs0[nc];
				}
			}
			
			double tv = 0, avdf = 0;
			if (svdf > 1.E-6) { tv = svt/svdf; avdf = svdf / sn;}

			cell_temp = (3.0 * tt + (srdf / sn) * trot + avdf * tv)/(3.0 + srdf / sn + avdf);
			cell_density = sn;
			if (tt < 0 ) {printf("!$! TT =  %e, sm*uu = %e\n", tt, (sm*uu/(3.0*K*sn)));}
			fprintf(plot_file, "%e %.6f; %.6f; %.6f %.6f; %.6f; %.6f; %.6f\n", num*dtp, tt, trot, tv, cell_temp, cell_density, alfa, denn);
			printf("tr = %.5f;  rot = %.5f; vib = %.5f; svib = %.5f; all = %.5f; alfa = %f; den = %e\n", tt, trot, tv, svib[0], cell_temp, alfa, den);
		}
	}
}

void show_temp_internal2(int num)
{
	int nc, x, y, krp, nm;
	double  norm;
	double tleft, tright, norml, normr;
	PARTICLE *p;
	VECTOR v0;
	nc = 0;
	for (y = 0; y < KY ; y++)
	{
		for (x = 0; x < KX; x++)
		{
			norm = 0.0;
			v0.x = v0.y = v0.z = 0.0;
			norml = tleft = 0.0;
			normr = tright = 0.0;
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						v0.x += p->v.x * p->weight;
						v0.y += p->v.y * p->weight;
						v0.z += p->v.z * p->weight;
						norm += p->weight;
					}
				}
				multiply(&v0, &v0, (1.0/norm));
			}
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						tleft += get_temp(p, &v0);
						tright += p->rot_energy*p->weight;
					}
				}
			}
			tleft /= norm;
			tright = tright * 4.0 / (norm * 2.0 * 2.0*K);
			fprintf(plot_file, "%.6f %.6f; %.6f; %.6f; %.6f\n", num*dtp,tleft, tright, (tleft * 3.0 + 2.0*tright) / (3.0 + 2.0), num*dtp*get_freq(tleft));
			printf("trans = %.5f;  rotational = %.5f; overal = %.5f\n", tleft, tright, (tleft * 3.0 + 2.0*tright) / (3.0 + 2.0));
		}
	}
}

void show_q_2flows(int num)
{
	FILE *fp;
	char file_name[60];
	int nc, x, y, krp, nm;
	double res_qflow, norm;
	double tleft, tright, norml, normr;
	PARTICLE *p;
	VECTOR v0;
	nc = 0;
	res_qflow = 0.0;
	sprintf(file_name, "%s/%i/%i%s", "temp2/q2flows", nc, num, ".txt");
	fp = fopen(file_name, "w+");
	for (y = 0; y < KY ; y++)
	{
		for (x = 0; x < KX; x++)
		{
			norm = 0.0;
			v0.x = v0.y = v0.z = 0.0;
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						v0.x += p->v.x * p->weight;
						v0.y += p->v.y * p->weight;
						v0.z += p->v.z * p->weight;
						norm += p->weight;
					}
				}
			}
			multiply(&v0, &v0, (1.0/norm));
			norml = tleft = 0.0;
			normr = tright = 0.0;
				
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					/* 					printf("v0.x = %.5f\n", v0.x); */
					for (krp = 0; krp < KRP; krp++)
					{
						p = &camera[krp][x][y][nc][nm];
						if (p->v.x >= 0)
						{
							tleft += (p->v.x - v0.x) * get_temp(p, &v0) * 1.5 * K;
							norml += p->weight;
						}
						if (p->v.x < 0)
						{
							tright += (p->v.x - v0.x) * get_temp(p, &v0) * 1.5 * K;
							normr += p->weight;
						}
					}
				}
			}
			res_qflow += (tleft+tright)/(norml+normr);
			fprintf(fp, "%.5f;   %.5f;		%.5f\n", (tleft / norml), (tright / normr), (tleft+tright)/(norml+normr));
		}
	}
	fclose(fp);
	printf("Q flow = %.5f\n", (res_qflow/KX));
}

void show_conc_2flows(int num)
{
	FILE *fp;
	char file_name[60];
	int nc, x, y, krp, nm;
	double cell_temp, norm;
	for (nc = 0; nc < NC; nc++)
	{
		sprintf(file_name, "%s/%i/%i%s", "temp2/conc2flows", nc, num, ".txt");
		fp = fopen(file_name, "w+");
		for (y = KY - 1; y >= 0 ; y--)
		{
			for (x = 0; x < KX; x++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					cell_temp = 0;
					norm = 0;
					for (krp = 0; krp < KRP; krp++)
					{
						cell_temp += camera[krp][x][y][nc][nm].weight;
					}
					fprintf(fp, "%.5f   ", (cell_temp/KRP));
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
		
	}
}

double get_cell_conc(int nc, int x, int y)
{
	int krp, nm;
	/*double fvel = 446.11, dvel;
	double uss = sqrt(5.0*K*temp[0]/(3.0*mass[0]));
	double ms = fvel/uss;
	double n1n2 = (ms*ms + 3.0)/(4.0*ms*ms);

	double n1 = conc[0]/KX;
	double n2 = n1n2*conc[0]/KX;*/
	double result = 0;
	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			result += camera[krp][x][y][nc][nm].weight;
			if (camera[krp][x][y][nc][nm].weight < 0)
			{
				printf("!!!   weight < 0   !!! (%i; %i)", x, y);
				printf("w=%e; v.x = %e; v.y = %e, v.z = %e; d = %e\n", camera[krp][x][y][nc][nm].weight, camera[krp][x][y][nc][nm].v.x,camera[krp][x][y][nc][nm].v.y,camera[krp][x][y][nc][nm].v.z, camera[krp][x][y][nc][nm].diam);
			}
		}
	}
	return (result/krp);
}

double get_cell_presure(int nc, int x, int y)
{
	int krp, nm;
	double result = 0, norm = 0;
	VECTOR v0;
	PARTICLE *p;
	v0.x = v0.y = v0.z = 0;
	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			p = &camera[krp][x][y][nc][nm];
			v0.x += p->v.x * p->weight;
			v0.y += p->v.y * p->weight;
			v0.z += p->v.z * p->weight;
			norm += p->weight;
		}
	}
	multiply(&v0, &v0, (1.0 / norm));

	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			result += get_temp(&camera[krp][x][y][nc][nm], &v0);
		}
	}
	return (K*result/krp);
}

void get_cell_impulse(VECTOR *v,int nc, int startx, int endx,  int starty, int endy)
{
	int krp, nm;
	int x, y;
	double norm = 0.0;
	v->x = 0;
	v->y = 0;
	v->z = 0;
	PARTICLE *p;
	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			for (x = startx; x <= endx; x++)
			{
				for (y = starty; y <= endy; y++)
				{
					p = &camera[krp][x][y][nc][nm];
					v->x += p->v.x*p->weight;
					v->y += p->v.y*p->weight;
					v->z += p->v.z*p->weight;
					norm +=	p->weight;
				}
			}
		}
	}
	multiply(v, v, 1.0/(norm));
}

void show_vector(int num)
{
	
	FILE *fp;
	char file_name[60];
	int nc, x, y;
	double cell_x, cell_y;
	VECTOR vec, sum;
	for (nc = 0; nc < NC; nc++)
	{
		sum.x = 0;
		sum.y = 0;
		sum.z = 0;
		sprintf(file_name, "%s/%i/%i%s",IMPULSE_DIR, nc, num, ".txt");
		fp = fopen(file_name, "w+");
		for (y = 0; y < KY ; y++)
		{
			for (x = 0; x < KX; x++)
			{
				get_cell_impulse(&vec, nc, x, x, y, y);
				cell_x = x * SIZEX;
				cell_y = y * SIZEY;
				fprintf(fp, "%.5f %.5f %.5f %.5f %.5f %.5f\n", cell_x, cell_y, vec.x * dtp, vec.y * dtp, vec.z * dtp, atan2(vec.y, vec.x));
				plus(&sum , &sum, &vec);
			}
		}
		printf("comp%i (%f; %f; %f)\n", nc, sum.x, sum.y, sum.z);
		fclose(fp);
	}
}

double get_cell_mass(int nc, int x, int y, int curr_mass)
{
	int krp, nm;
	PARTICLE *p;
	double result = 0;
	for(krp=0; krp<KRP; krp++)
	{
		for (nm=0; nm<NM; nm++)
		{
			p = &camera[krp][x][y][nc][nm];
			if ((p->mass>=curr_mass*ETALON_MASS)&&(p->mass<(curr_mass+1)*ETALON_MASS))
			{
				result += p->weight;
			}
		}
	}
	return (result/(KRP));
}

double get_cell_full_mass(int nc, int x, int y)
{
	int krp, nm;
	PARTICLE *p;
	double result = 0, summ = 0;
	for(krp=0; krp<KRP; krp++)
	{
		for (nm=0; nm<NM; nm++)
		{
			p = &camera[krp][x][y][nc][nm];
			result += p->weight*p->mass;
			summ += p->weight;
		}
	}
	if (0 == summ) return 0;
	return (result/summ);
}

void show_mass(int num)
{
	FILE *fp;
	char file_name[60];
	int nc, x, y, curr_mass;
	double full_mass[KX][KY], sum[PIC_MASS_COUNT], overal_mass;
	overal_mass = 0;
	fprintf(plot_file, "\n%e; ", num * dtp);
	for (nc = 0; nc < NC; nc++)
	{
		for (curr_mass = 0; curr_mass < PIC_MASS_COUNT; curr_mass++)
		{
			
			for (y=0; y<KY; y++)
			{
				for (x=0; x<KX; x++)
				{
					full_mass[x][y] = get_cell_mass(nc, x, y, (int)pow(10,curr_mass));
				}
			}
			sum[curr_mass] = 0;
			sprintf(file_name, "%s/%i/%i%s",MASS_DIR, curr_mass, num, ".txt");
// 			printf("\n%s\nOk\n",file_name);
			fp = fopen(file_name, "w+");
			for (y = 0; y < KY; y++)
			{
				for (x = 0; x < KX; x++)
				{
					sum[curr_mass] += full_mass[x][y];
					fprintf(fp, "%.5f ", (full_mass[x][y]));
				}
				fprintf(fp, "\n");
			}
 			overal_mass += sum[curr_mass] * curr_mass * mass[nc];
			printf("mass%i = %e, avr = %e; ", (int)pow(10,curr_mass), sum[curr_mass]/(conc[nc]), sum[curr_mass] / (KX * KY * conc[nc]));
			fprintf(plot_file,"%e; ", sum[curr_mass] / (KX * KY * conc[nc]));
			fclose(fp);
		}
		printf("\nmass = %e;\n", overal_mass);
	}
}

void show_nalip(int num)
{
	FILE *fpl, *fpr, *fpu, *fpd;
	char file_left[60], file_right[60], file_down[60], file_up[60];
	int nc, x, y;
	double full_left[KY][NC][NALIP_MASS_COUNT], full_right[KY][NC][NALIP_MASS_COUNT];
	double full_up[KX][NC][NALIP_MASS_COUNT], full_down[KX][NC][NALIP_MASS_COUNT];
	for (nc = 0; nc < NC; nc++)
	{
		
		sprintf(file_left, "%s/%i/%s/%i%s", NALIP_DIR, nc, "left", num, ".txt");
		sprintf(file_right, "%s/%i/%s/%i%s", NALIP_DIR, nc, "right", num, ".txt");
		sprintf(file_up, "%s/%i/%s/%i%s", NALIP_DIR, nc, "up", num, ".txt");
		sprintf(file_down, "%s/%i/%s/%i%s", NALIP_DIR, nc, "down", num, ".txt");
		/*printf("%s\n", file_name);*/
		fpl = fopen(file_left, "w+");
		fpr = fopen(file_right, "w+");
		fpd = fopen(file_down, "w+");
		fpu = fopen(file_up, "w+");
		for (x=0; x<NALIP_MASS_COUNT; x++)
		{
			for (y=0; y<KY; y++)
			{
				fprintf(fpl, "%.5f ", (full_left[y][nc][x]/((num_proc-1)*KRP)));
				fprintf(fpr, "%.5f ", (full_right[y][nc][x]/((num_proc-1)*KRP)));
			}
			fprintf(fpl, "\n");
			fprintf(fpr, "\n");
		}
		fclose(fpl);
		fclose(fpr);
		for (x=0; x<NALIP_MASS_COUNT; x++)
		{
			for (y=0; y<KX; y++)
			{
				fprintf(fpu, "%.5f ", (full_up[y][nc][x]/((num_proc-1)*KRP)));
				fprintf(fpd, "%.5f ", (full_down[y][nc][x]/((num_proc-1)*KRP)));
			}
			fprintf(fpu, "\n");
			fprintf(fpd, "\n");
		}
		fclose(fpu);
		fclose(fpd);
	}
}

void show_reaction_rates(int num)
{
	int krp, i, j;
	double kab = 0.0, kcd = 0.0;
	PARTICLE *pl, *pm;
	VECTOR vrc;
	double g, sigma, min_weight;
	double kab_eq=0.0, kcd_eq=0.0;
	double t_ab = 0.0, t_cd = 0.0;
	double m_ab = 0.5, m_cd = 0.5;
	double n_a, n_b, n_c, n_d;
	double kab_ij=0.0, kcd_ij=0.0;
	n_a = get_cell_conc(0,0,0);
	n_b = get_cell_conc(1,0,0);
	n_c = get_cell_conc(2,0,0);
	n_d = get_cell_conc(3,0,0);
	t_ab = (get_t(0,0,0)+get_t(1,0,0))/2.0;
	t_cd = (get_t(2,0,0)+get_t(3,0,0))/2.0;
	for (krp = 0; krp < KRP; krp++)
	{
		for (i = 0; i < NM; i++)
		{
			for (j = 0; j < NM; j++)
			{
				pl = &camera[krp][0][0][0][i];
				pm = &camera[krp][0][0][1][j];
				min_weight = (pl->weight<pm->weight? pl->weight :pm->weight);
				
				minus(&vrc, &(pl->v), &(pm->v));
				g = mag(&vrc);
 				sigma = PI*pow(pl->diam+pm->diam, 2)*(1.0 - 4.0*3.0*K/(g*g))/4.0;

				if (sigma>0)
				{
					kab_ij += sigma*g*pl->weight*pm->weight;
				}
				
				pl = &camera[krp][0][0][2][i];
				pm = &camera[krp][0][0][3][j];
				min_weight = (pl->weight<pm->weight? pl->weight :pm->weight);
				minus(&vrc, &(pl->v), &(pm->v));
				g = mag(&vrc);
				sigma = PI*pow(pl->diam+pm->diam, 2)*(1.0 - 4.0*53.0*K/(g*g))/4.0;

				if (sigma>0)
				{
					kcd_ij += sigma*g*pl->weight*pm->weight;
				}
			}
		}
	}
	sigma = pow(camera[0][0][0][0][0].diam, 2);
	kab_eq += 2.0 *sigma* sqrt(2.0*PI*K*t_ab/m_ab)*exp(-3.0/(t_ab));
	kcd_eq += 2.0 *sigma* sqrt(2.0*PI*K*t_cd/m_cd)*exp(-53.0/(t_cd));
	
	kab = kab_ij/(n_a*n_b*KRP);
	kcd = kcd_ij/(n_c*n_d*KRP);
/* 	kab = kab_ij / (na); */
/* 	kcd = kcd_ij / (nc);
	double g, sigma, min_weight; */

	printf("Kab = %8.8f; Kcd = %8.8f; Kab_eq = %8.8f; Kcd_eq = %8.8f; Kf = %8.8f; Kb = %8.8f\n", kab, kcd, kab_eq, kcd_eq, kab/kab_eq, kcd/kcd_eq);
	
	FILE *fp;
	char file_name[60];
	sprintf(file_name, "temp2/rates/%i%s", num, ".txt");
	fp = fopen(file_name, "w+");
	fprintf(fp, "%8.8f %8.8f", kab/kab_eq, kcd/kcd_eq);
	fclose(fp);
}

void show_dissoc_rate(int step)
{
	int krp, nm;
	double n_a;
	double kab_ij=0.0, k_av = 0.0;

	n_a = 0;
	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			n_a += camera[krp][0][0][0][nm].weight;
		}
	}

	total_dissoc += step_dissoc/KRP;
	
	//t_ab = get_t(0,0,0);
	
	if (n_a>0) {
		kab_ij = (step_dissoc / dtp)/(n_a*n_a);
		k_av = (total_dissoc / (step * dtp))/(n_a*n_a);
	}

	
	//kab_eq += 2.0 *sigma* sqrt(2.0*PI*K*t_ab/m_ab)*exp(-3.0/(t_ab));
	//kcd_eq += 2.0 *sigma* sqrt(2.0*PI*K*t_cd/m_cd)*exp(-53.0/(t_cd));
	
	printf("Kav = %e; Kmg = %e;  dissoc = %e; total_dissoc = %e\n", k_av, kab_ij, (step_dissoc/KRP), total_dissoc);
	step_dissoc = 0;
}

void save_field(char* file_dir, double*** field_array, int time_step)
{
	int kx, ky, nc;
	FILE *fp;
	char file_name[60];
	for (nc = 0; nc  < NC; nc++)
	{
		sprintf(file_name, "%s/%i/%i%s", file_dir, nc, time_step, ".txt");
		fp = fopen(file_name, "w+");
		for (ky = 0; ky < KY; ky++)
		{
			for (kx = 0; kx < KX; kx++)
			{
				fprintf(fp, "%.8e; ", field_array[nc][kx][ky]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
}

void save_three_column(char* file_dir, double*** field_array, int time_step)
{
	int kx, ky, nc;
	FILE *fp;
	char file_name[60];
	for (nc = 0; nc  < NC; nc++)
	{
		sprintf(file_name, "%s/%i/%i%s", file_dir, nc, time_step, ".txt");
		fp = fopen(file_name, "w+");
		for (ky = 0; ky < KY; ky++)
		{
			for (kx = 0; kx < KX; kx++)
			{
				fprintf(fp, "%4.4e; %4.4e; %4.4e; %4.4e\n", kx*SIZEX, ky*SIZEY, field_array[nc][kx][ky], temp_field[nc][kx][ky]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
}

void save_plot(char* file_dir, double** field_array, int time_step)
{
	int nc, i;
	FILE *fp;
	char file_name[60];
	for (nc = 0; nc  < NC; nc++)
	{
		sprintf(file_name, "%s/%i/%i%s", file_dir, nc, time_step, ".txt");
		fp = fopen(file_name, "w+");
		for (i = 0; i < MAX_PLOT_VEL; i++)
		{
			fprintf(fp, "%.8e; %.8e;\n", i - 0.5*MAX_PLOT_VEL,field_array[nc][i]);
		}
		fclose(fp);
	}
}

void print_forces(int time_step)
{
	int x, y; double mag_f;
	FILE *fx, *fy, *fm;
	char file_name_x[60], file_name_y[60], file_name_m[60];
	sprintf(file_name_x, "temp2/force/x/%i.txt", time_step);
	fx = fopen(file_name_x, "w+");
	sprintf(file_name_y, "temp2/force/y/%i.txt", time_step);
	fy = fopen(file_name_y, "w+");
	sprintf(file_name_m, "temp2/force/m/%i.txt", time_step);
	fm = fopen(file_name_m, "w+");
	for (y = 0; y < KY; y++)
	{
		for (x = 0; x < KX; x++)
		{
			fprintf(fx, "%.8e; ", force_field[x][y].x/KRP);
			fprintf(fy, "%.8e; ", force_field[x][y].y/KRP);
			mag_f = sqrt(pow(force_field[x][y].x/KRP,2)+pow(force_field[x][y].y/KRP,2));
			fprintf(fm, "%.8e; ", mag_f);
			force_field[x][y].x = 0;
			force_field[x][y].y = 0;
			force_field[x][y].z = 0;
		}
		fprintf(fx, "\n");
		fprintf(fy, "\n");
		fprintf(fm, "\n");
	}
	fclose(fx);
	fclose(fy);
	fclose(fm);
}

void bird_style_pics(int step)
{
	int kx, ky, krp, nc, nm;
	PARTICLE *p;
	double avr_conc[NC], avr_temp[NC], avr_mass[NC], new_avr_vx[NC], new_avr_vy[NC];
	for (nc = 0; nc < NC; nc++)
	{
		new_avr_vx[nc] = 0;
		new_avr_vy[nc] = 0;
		for (kx = 0; kx < KX; kx++)
		{
			for (ky = 0; ky < KY; ky++)
			{
				temp_field[nc][kx][ky] = 0.0;
				conc_field[nc][kx][ky] = 0.0;
				mass_field[nc][kx][ky] = 0.0;
				impulse_field[nc][kx][ky] = 0.0;
			}
		}
		for (kx = 0; kx < MAX_PLOT_VEL; kx++)
		{
			distr_x[nc][kx] = 0;
			distr_y[nc][kx] = 0;
		}
	}
	for (krp = 0; krp < KRP; krp++)
	{
		for (kx = 0; kx < KX; kx++)
		{
			for (ky = 0; ky < KY; ky++)
			{
				for (nc = 0; nc < NC; nc++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						p = &camera[krp][kx][ky][nc][nm];
						conc_field[nc][kx][ky] += p->weight;
						temp_field[nc][kx][ky] += get_energy(p);
						mass_field[nc][kx][ky] += p->weight * p->mass * SIZEX * SIZEY * SIZEZ;
						impulse_field[nc][kx][ky] += p->weight * atan2(p->v.y, p->v.x);
						new_avr_vx[nc] += p->weight * p->v.x;
						new_avr_vy[nc] += p->weight * p->v.y;
						int num_vel_x = (int)(p->v.x + 0.5 * MAX_PLOT_VEL), num_vel_y = (int)(p->v.y + 0.5 * MAX_PLOT_VEL);
						if ((num_vel_x < MAX_PLOT_VEL)&&(num_vel_x>0)) {distr_x[nc][num_vel_x] += p->weight;}
						if ((num_vel_y < MAX_PLOT_VEL)&&(num_vel_y>0)) {distr_y[nc][num_vel_y] += p->weight;}
// 						if (fabs(p->v.x) < 0.5*MAX_PLOT_VEL) {distr_x[nc][(int)(p->v.x + 0.5 * MAX_PLOT_VEL)] += p->weight;}
// 						if (fabs(p->v.y) < 0.5*MAX_PLOT_VEL) {distr_y[nc][(int)(p->v.y + 0.5 * MAX_PLOT_VEL)] += p->weight;}
					}
				}
			}
		}
	}
	for (nc = 0; nc < NC; nc++)
	{
		avr_conc[nc] = 0;
		avr_temp[nc] = 0;
		avr_mass[nc] = 0;
		for (kx = 0; kx < KX; kx++)
		{
			for (ky = 0; ky < KY; ky++)
			{
				temp_field[nc][kx][ky] /= KRP;
				avr_temp[nc] += temp_field[nc][kx][ky];
				conc_field[nc][kx][ky] /= KRP;
				avr_conc[nc] += conc_field[nc][kx][ky];
				mass_field[nc][kx][ky] /= KRP;
				avr_mass[nc] += mass_field[nc][kx][ky];
				impulse_field[nc][kx][ky] = (conc_field[nc][kx][ky]>0)?(impulse_field[nc][kx][ky]/(conc_field[nc][kx][ky]*KRP)):100.0;
			}
		}
		avr_vx[nc] = (avr_conc[nc] > 0) ? new_avr_vx[nc]/(avr_conc[nc] * KRP) : 0;
		avr_vy[nc] = (avr_conc[nc] > 0) ? new_avr_vy[nc]/(avr_conc[nc] * KRP) : 0;
		printf("comp#%i: temp = %5.5e; conc = %5.5e; mass = %5.5e; vx = %e; vy = %e;\n", nc, avr_temp[nc]/(KX*KY), avr_conc[nc]/(KX*KY), avr_mass[nc] /(KX*KY), avr_vx[nc], avr_vy[nc]);
		
	}
	save_plot("temp2/plotx", distr_x, step);
	save_plot("temp2/ploty", distr_y, step);
	save_field("temp2/temp", temp_field, step);
	save_field("temp2/conc", conc_field, step);
	save_field("temp2/full_mass", mass_field, step);
	save_three_column("temp2/impulse", impulse_field, step);
//	print_forces(step);
}

void print_pictures(int step)
{
	if ((0 != SHOW_TEMP_FREQ) && !(step % SHOW_TEMP_FREQ)) bird_style_pics(step);
//	printf("fin");
//   	if ((0 != SHOW_TEMP_FREQ) && !(step % SHOW_TEMP_FREQ)) show_field(step, get_t, "temp2/temp", 0, KX-1, 0, KY-1);
//   	fprintf(plot_file,"%e; %e; %e; \n", step * dtp, get_t(0,0,0), get_t(1,0,0));

	/* 	if ((0 != SHOW_CONC_FREQ) && !(step % SHOW_CONC_FREQ)) show_field(step, get_qflow, "temp2/qflow"); */

//  	if ((0 != SHOW_TEMP_FREQ) && !(step % SHOW_TEMP_FREQ)) show_field(step, get_cell_full_mass, "temp2/full_mass", 0, KX-1, 0, KY-1);

 	//if ((0 != SHOW_IMPULSE_FREQ) && !(step % SHOW_IMPULSE_FREQ)) show_vector(step);

}

