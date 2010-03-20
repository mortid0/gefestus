#include "params.h"
#include <math.h>
#include "utils.h"
#include <stdio.h>
#define NUM_SAMP 9

void normirovka() 
{
	double c[NC][6],sn = 0, snt = 0, a = 0, tfe;
	int i,j;

	for (i = 0; i < NC; i++)
	{
		printf("diam[%i] = %f\n", i, diam[i]);
	}
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < 6; j++)
		{
			c[i][j] = 0;
		}
	}
	for (i = 0; i < NC; i++)
	{
		c[i][0] = conc[i];
		c[i][3] = conc[i]*temp[i];
		c[i][1] = mass[i];
		c[i][2] = diam[i];
	}
	if (IS_INJECT)
	{
		for (i = 0; i < NC; i++)
		{
			c[i][0] += CONC_INJ;
			c[i][3] += CONC_INJ * TEMP_INJ;
		}
	}
	for (i = 0; i < NC; i++)
	{
		sn += c[i][0];
		snt += c[i][0]*c[i][3];
		c[i][2] = c[i][2]/c[0][2];
		c[i][4] = c[i][0];
		c[i][5] = 0;
	}
	tfe = snt / sn;
	// 	tfe = 1.0;
	printf("tfe = %.5f\n", tfe);
	for (i = 0; i < NC; i++)
	{
		for (j = 0; j < NC; j++)
		{
			a += (c[i][0]/sn)*(sqrt(PI)/2)*(pow(c[i][2]+c[j][2],2)*c[j][0]*sqrt(tfe*(c[i][1]+c[j][1])/(c[i][1]*c[j][1])));
		}
	}
	c[0][2] = 1.0/sqrt(a);
	diam[0] = c[0][2];
	for (i = 1; i < NC; i++)
	{
		c[i][2] = c[0][2]*c[i][2];
		diam[i] = c[i][2];
	}
	for (i = 0; i < NC; i++)
	{
		printf("diam[%i] = %f\n", i, diam[i]);
	}
}

void mixture()
{
	int krp, x, y, nc, nm, r_krp, r_nm;
	PARTICLE temp_particle, *p1, *p2;

	for (nc = 0; nc < NC; nc++)
	{
		for (x = 0; x < KX; x++)
		{
			for (y = 0; y < KY; y++)
			{
				for (krp = 0; krp < KRP; krp++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						r_krp = rand() % KRP;
						r_nm = rand() % NM;
						/*if (r_krp >= KRP) printf("krp = %i\n", r_krp);
						  if (r_nm >= NM) printf("nm = %i\n", r_nm);*/
						p1 = &camera[krp][x][y][nc][nm];
						p2 = &camera[r_krp][x][y][nc][r_nm];
						copy_particle(p1, &temp_particle);
						copy_particle(p2, p1);
						copy_particle(&temp_particle, p2);
					}
				}
			}
		}
	}
}

void inject(double time)
{
	double d, a;
	int i, nc, krp;
	if (-1 == CELL_INJ) return;
	if (time > TIME_INJ) return;
	d = sqrt(TEMP_INJ/MASS_INJ);
	for (krp = 0; krp < KRP; krp++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			for (i = 0; i < NM; i++)
			{
				camera[krp][0][CELL_INJ][nc][i].weight += CONC_INJ * dtp /(NM * TIME_INJ);
				a = CONC_INJ / (NM * camera[krp][0][CELL_INJ][nc][i].weight);
				if (frand() < a) 
				{
					camera[krp][0][CELL_INJ][nc][i].v.x = d*sqrt(-log(frand()));
					camera[krp][0][CELL_INJ][nc][i].v.y = d*sqrt(-log(frand()))*cos(2*PI*frand());
					camera[krp][0][CELL_INJ][nc][i].v.z = d*sqrt(-log(frand()))*sin(2*PI*frand());
					camera[krp][0][CELL_INJ][nc][i].mass = MASS_INJ;
					camera[krp][0][CELL_INJ][nc][i].diam = DIAM_INJ;
				}
			} 
		}
	}
}

void inject_left_right_half(double time)
{
	double d, a, vr, phi;
	int i, nc, krp, y;

	if (time > TIME_INJ) return;
	d = sqrt(TEMP_INJ/MASS_INJ);

	for (krp = 0; krp < KRP; krp++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			for (y = 0; y < KY/2; y++)
			{
				for (i = 0; i < NM; i++)
				{
					camera[krp][0][y][nc][i].weight += CONC_INJ * dtp /(NM * TIME_INJ);
					a = CONC_INJ / (NM * camera[krp][0][y][nc][i].weight);
					if (frand() < a) 
					{
						camera[krp][0][y][nc][i].v.x = d;
						vr = d*sqrt(-log(frand()));
						phi = 2.0*PI*frand();
						camera[krp][0][y][nc][i].v.y = vr*cos(phi);
						camera[krp][0][y][nc][i].v.z = vr*sin(phi);
						camera[krp][0][y][nc][i].mass = MASS_INJ;
						camera[krp][0][y][nc][i].diam = DIAM_INJ;
					}
				} 
			}
		}
	}
	for (krp = 0; krp < KRP; krp++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			for (y=KY/2; y< KY; y++)
			{
				for (i = 0; i < NM; i++)
				{
					camera[krp][KX-1][y][nc][i].weight += CONC_INJ * dtp /(NM * TIME_INJ);
					a = CONC_INJ / (NM * camera[krp][KX-1][y][nc][i].weight);
					if (frand() < a) 
					{
						camera[krp][KX-1][y][nc][i].v.x = -d;
						vr = d*sqrt(-log(frand()));
						phi = 2.0*PI*frand();
						camera[krp][KX-1][y][nc][i].v.y = vr*cos(phi);
						camera[krp][KX-1][y][nc][i].v.z = vr*sin(phi);

						camera[krp][KX-1][y][nc][i].mass = MASS_INJ;
						camera[krp][KX-1][y][nc][i].diam = DIAM_INJ;
					}
				} 
			}
		}
	}
}

void inject_left(double time)
{
	double d, a;
	int i, nc, krp, y;
	double vr, phi;
	if (time > TIME_INJ) return;
	d = sqrt(2.0*K*TEMP_INJ/MASS_INJ);

	for (krp = 0; krp < KRP; krp++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			for (y = (KY/2); y < KY; y++)
			{
				for (i = 0; i < NM; i++)
				{
					camera[krp][0][y][nc][i].weight += CONC_INJ * dtp /(NM * TIME_INJ);
					a = CONC_INJ / (NM * camera[krp][0][y][nc][i].weight);
					if (frand() < a) 
					{
						camera[krp][0][y][nc][i].v.x = d*sqrt(-log(frand()));
						vr = d*sqrt(-log(frand()));
						phi = 2.0*PI*frand();
						camera[krp][0][y][nc][i].v.y = vr*sin(phi);
						camera[krp][0][y][nc][i].v.z = vr*cos(phi);
						camera[krp][0][y][nc][i].mass = MASS_INJ;
						camera[krp][0][y][nc][i].diam = DIAM_INJ;
					}
				} 
			}
		}
	}
}

void inject_right(double time)
{
	double d, a;
	int i, nc, krp, y;
	double vr, phi;
	if (time > TIME_INJ) return;
	d = sqrt(2.0*K*TEMP_INJ/MASS_INJ);

	for (krp = 0; krp < KRP; krp++)
	{
		for (nc = 0; nc < NC; nc++)
		{
			for (y = 0; y < KY; y++)
			{
				for (i = 0; i < NM; i++)
				{
					camera[krp][KX-1][y][nc][i].weight += CONC_INJ * dtp /(NM * TIME_INJ);
					a = CONC_INJ / (NM * camera[krp][KX-1][y][nc][i].weight);
					if (frand() < a) 
					{
						camera[krp][KX-1][y][nc][i].v.x = -d*sqrt(-log(frand()));
						vr = d*sqrt(-log(frand()));
						phi = 2.0*PI*frand();
						camera[krp][KX-1][y][nc][i].v.y = vr*sin(phi);
						camera[krp][KX-1][y][nc][i].v.z = vr*cos(phi);
						camera[krp][KX-1][y][nc][i].mass = MASS_INJ;
						camera[krp][KX-1][y][nc][i].diam = DIAM_INJ;
					}
				} 
			}
		}
	}
}

void change_plate(int krp)
{

}

void accurate_force(int krp)
{
	int x, y, i, j, nm, nc;
	double m[KX][KY], rad, f12, phi, vol = SIZEX*SIZEY;
	double f_x, f_y;
	VECTOR accl[KX][KY];
	PARTICLE *p;
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			m[x][y] = 0;
			accl[x][y].x = 0;
			accl[x][y].y = 0;
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					m[x][y] += camera[krp][x][y][nc][nm].mass*camera[krp][x][y][nc][nm].weight;
				}
			}
		}
	}
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			for (i = x; i < KX; i++)
			{
				for (j = y; j < KY; j++)
				{
					if ((i ==x ) && (j == y)) {continue;}
					rad = sqrt((pow((j-y)*SIZEY, 2) + pow((i-x)*SIZEX, 2)));
					f12 = GAM * vol / (rad);
					phi = atan2((y-j)*SIZEY,((i - x)*SIZEX));
					f_x = f12 * cos(phi);
					f_y = f12 * sin(phi);
					accl[x][y].x += f_x*m[i][j];
					accl[x][y].y += f_y*m[i][j];

					accl[i][j].x -= f_x*m[x][y];
					accl[i][j].y -= f_y*m[x][y];
				}
			}
		}
		// 	printf("Calculation force finished\n");
		for (x = 0; x < KX; x++)
		{
			for (y = 0; y < KY; y++)
			{
				for (nc = 0; nc < NC; nc++)
				{
					for (nm = 0; nm < NM; nm++)
					{
						p = &camera[krp][x][y][nc][nm];
						if (0 == p->weight) {continue;}
						p->v.x += accl[x][y].x*dtp;
						p->v.y += accl[x][y].y*dtp;
					}
				}
			}
		}

	}
}

void force(int krp)
{
	int x, y, i, j, nm, nc, num_samples;
	double m[KX][KY], cum_m[KX*KY], rad, f12, phi, cam_mass = 0, vol = SIZEX*SIZEY;
	VECTOR accl[KX][KY];
	PARTICLE *p;
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			accl[x][y].x = 0;
			accl[x][y].y = 0;
			m[x][y] = 0;
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					m[x][y] += camera[krp][x][y][nc][nm].mass*camera[krp][x][y][nc][nm].weight;
				}
			}
			if ((0==x) && (0 ==y)) {cum_m[0] = m[x][y];}
			else {cum_m[x*KY+y] = cum_m[x*KY+y-1] + m[x][y];}
		}
	}
	cam_mass = cum_m[KX*KY-1];
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			num_samples = 0;
			double ver;int found = 0, sred;
			while (num_samples < NUM_SAMP){
				ver = frand() * (cam_mass - m[x][y]);
				if (ver > cum_m[x*KY+y - 1]) {ver+= m[x][y];}
				i = 0;
				j = KY*KX - 1;
				found = -1;
				if (ver < cum_m[0]) {found = 0;}
				if (ver > cum_m[KX*KY-2]) {found = j;}
				do {
					sred = (i-j)/2 + j;
					if ((cum_m[sred]>ver)&&(cum_m[sred-1]<ver)) {found = sred; break;}
					else if (cum_m[sred]<ver) {i = sred;}
					else {j = sred;}

				} while ((j>i)&&(-1==found));
				i = found / KY; j = found % KY;
				num_samples++;

				rad = sqrt((pow((j-y)*SIZEY, 2) + pow((i-x)*SIZEX, 2)));
				f12 = GAM * vol / (rad);
				phi = atan2((y-j)*SIZEY,((i - x)*SIZEX));

				accl[x][y].x += f12 * (cam_mass - m[x][y]) * cos(phi) / NUM_SAMP;
				accl[x][y].y += f12 * (cam_mass - m[x][y]) * sin(phi) / NUM_SAMP;
				accl[i][j].x -= f12 * m[x][y] * cos(phi);
				accl[i][j].y -= f12 * m[x][y] * sin(phi);
			}
		}
	}
	// 	printf("Calculation force finished\n");
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			for (nc = 0; nc < NC; nc++)
			{
				for (nm = 0; nm < NM; nm++)
				{
					p = &camera[krp][x][y][nc][nm];
					if (0 == p->weight) {continue;}
					p->v.x += accl[x][y].x*dtp;
					p->v.y += accl[x][y].y*dtp;
				}
			}
		}
	}
}
