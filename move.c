#include "move.h"
#include "params.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#define Procrustes 1.0

void null_buf()
{
	int x, y, nm;
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			for (nm = 0; nm < NM; nm++)
			{
				buf[x][y][nm].weight = 0;
			}
		}
	}
}

void no_vertical_wall_move(int krp, int x, int y, int nc, int nm)
{
	double vx, vy, wt;
	double part1, part3;
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	int zy = ((p->v.y) > 0)?1:-1;
	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	p->weight = part1;
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &(p->v), part3);
}

void final_buf(int nc, int krp)
{
	int x, y, nm;
	PARTICLE *p ,*b;
// 	double rn;
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			for (nm = 0; nm < NM; nm++)
			{
				p = &(camera[krp][x][y][nc][nm]);
				b = &(buf[x][y][nm]);
				p->weight += (b->weight);
				if (frand() < (b->weight / p->weight))
				{
					p->v.x = b->v.x;
					p->v.y = b->v.y;
					p->v.z = b->v.z;
					p->diam = b->diam;
					p->mass = b->mass;
					p->rot_energy = b->rot_energy;
					p->rot_degr = b->rot_degr;
					//p->reflect_field = b->reflect_field;
				}
				//��������� �� �������. ���� �� ���� ��� �� ������� ������� ����� ��������� ������, ��� �������������, �� ��������� ��� �� ������� ��� �������� � 2 ����.
				if (p->v.x * dtp_var[krp] > SIZEX * kProlet)
				{
					dtp_var[krp] /= 4.0;
					printf("!!! Prolet [%i,%i], vx = %e; vy = %e; m = %e; d = %e; w = %e\n",x, y, p->v.x, p->v.y, p->mass, p->diam, p->weight);
					printf("new dtp_var[%i] = %f", krp, dtp_var[krp]);
				}
				if (p->v.y * dtp_var[krp] > SIZEY * kProlet)
				{
					dtp_var[krp] /= 4.0;
					printf("!!! Prolet [%i,%i], vx = %e; vy = %e; m = %e; d = %e; w = %e\n",x, y, p->v.x, p->v.y, p->mass, p->diam, p->weight);
					printf("new dtp_var[%i] = %f", krp, dtp_var[krp]);
				}
			}
		}
	}
}

void inside_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	double part1, part2, part3, part4;
	
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	zx = ((p->v.x) > 0)?1:-1;
	zy = ((p->v.y) > 0)?1:-1;

	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part2 = vx * dtp_var[krp] * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	part4 = vx * dtp_var[krp] * vy * dtp_var[krp] * wt;
	p->weight = part1;
	
	move_particle_to_buf(p, &(buf[x+zx][y][nm]), &(p->v), part2);
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &(p->v), part3);
	move_particle_to_buf(p, &(buf[x+zx][y+zy][nm]), &(p->v), part4);
}

void vertical_wall_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	VECTOR wallx;	
	double part1, part2, part3, part4;
	
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	zx = 0;
	zy = ((p->v.y) > 0)?1:-1;

	copy_vector(&p->v, &wallx);
	wallx.x *= (-1.0);
	
	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part2 = vx * dtp_var[krp] * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	part4 = vx * dtp_var[krp] * vy * dtp_var[krp] * wt;
	p->weight = part1;
	//p->reflect_field = 0;
	
	move_particle_to_buf(p, &(buf[x+zx][y][nm]), &wallx, part2);
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &(p->v), part3);
	move_particle_to_buf(p, &(buf[x+zx][y+zy][nm]), &wallx, part4);
}

void horizontall_wall_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	VECTOR wally;	
	double part1, part2, part3, part4;
	
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	zx = ((p->v.x) > 0)?1:-1;
	zy = 0;

	copy_vector(&p->v, &wally);
 	wally.y *= (-1.0);
// 	wally.x *= Procrustes;
// 	wally.z *= Procrustes;
	
	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part2 = vx * dtp_var[krp] * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	part4 = vx * dtp_var[krp] * vy * dtp_var[krp] * wt;
	p->weight = part1;
	//p->reflect_field = 0;
	move_particle_to_buf(p, &(buf[x+zx][y][nm]), &(p->v), part2);
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &wally, part3);
	move_particle_to_buf(p, &(buf[x+zx][y+zy][nm]), &wally, part4);
}

void no_horizontall_wall_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	VECTOR wally;	
	double part1, part2;
	
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	zx = ((p->v.x) > 0)?1:-1;
	zy = 0;

	copy_vector(&p->v, &wally);
	wally.y *= (-1.0);
	
	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part2 = vx * dtp_var[krp] * (SIZEY - vy * dtp_var[krp]) * wt;
	
	p->weight = part1;
	
	move_particle_to_buf(p, &(buf[x+zx][y][nm]), &(p->v), part2);	
}

void inside_angle_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	double part1, part2, part3, part4;
	
	PARTICLE *p;
	VECTOR wallx, wally, wallxy;
	p = &(camera[krp][x][y][nc][nm]);
	zx = 0;
	zy = 0;
	copy_vector(&(p->v), &wallx);
	wallx.x *= (-1.0);
	copy_vector(&(p->v), &wally);
	wally.y *= (-1.0);
	copy_vector(&wallx, &wallxy);
	wallxy.y *= (-1.0);

	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part2 = vx * dtp_var[krp] * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	part4 = vx * dtp_var[krp] * vy * dtp_var[krp] * wt;
	p->weight = part1;
	
	move_particle_to_buf(p, &(buf[x+zx][y][nm]), &wallx, part2);
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &wally, part3);
	move_particle_to_buf(p, &(buf[x+zx][y+zy][nm]), &wallxy, part4);
}

void inside_angle_no_vert_move(int krp, int x, int y, int nc, int nm)
{
	int zx,zy;
	double vx, vy, wt;
	double part1, part3;
	
	PARTICLE *p;
	VECTOR wallx, wally, wallxy;
	p = &(camera[krp][x][y][nc][nm]);
	zx = (p->v.x>0)?(1):(-1);
	zy = 0;
// 	if ((zx + x >= 0) || (zx + x < KX)){horizontall_wall_move(krp, x, y, nc, nm); return;}
	copy_vector(&(p->v), &wallx);
	wallx.x *= (-1.0);
	copy_vector(&(p->v), &wally);
	wally.y *= (-1.0);
	copy_vector(&wallx, &wallxy);
	wallxy.y *= (-1.0);

	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	part3 = (SIZEX - vx * dtp_var[krp]) * vy * dtp_var[krp] * wt;
	p->weight = part1;
	
	move_particle_to_buf(p, &(buf[x][y+zy][nm]), &wally, part3);
}

void inside_angle_no_wall_move(int krp, int x, int y, int nc, int nm)
{
	int zx, zy;
	double vx, vy, wt;
	double part1;
	
	PARTICLE *p;
	VECTOR wallx, wally, wallxy;
	p = &(camera[krp][x][y][nc][nm]);
	zx = 0;
	zy = 0;
	copy_vector(&(p->v), &wallx);
	wallx.x *= (-1.0);
	copy_vector(&(p->v), &wally);
	wally.y *= (-1.0);
	copy_vector(&wallx, &wallxy);
	wallxy.y *= (-1.0);
// 	multiply(&wally, &wally, Procrustes);

	vx = fabs(p->v.x);
	vy = fabs(p->v.y);
	wt = p->weight / (SIZEX * SIZEY);
	part1 = (SIZEX - vx * dtp_var[krp]) * (SIZEY - vy * dtp_var[krp]) * wt;
	p->weight = part1;
}

void make_2d_move(int krp, int x, int y, int nc, int nm)
{
	PARTICLE *p;
	p = &(camera[krp][x][y][nc][nm]);
	if (0 == p->weight) return;
	
	int zx = ((p->v.x) > 0)?1:-1;
	int zy = ((p->v.y) > 0)?1:-1;

	if (( KX == (x + zx)) || (-1 == (x + zx)))
	{
		if (KY == (y + zy) || (-1 == (y + zy)))
		{
 			inside_angle_no_wall_move(krp, x, y, nc, nm);
			return;
		}
     		no_vertical_wall_move(krp, x, y, nc, nm );
		return;
	}
	if (KY == (y + zy) || (-1 == (y + zy)))
	{
    		no_horizontall_wall_move(krp, x, y, nc, nm);
		return;
	}
 	inside_move(krp, x, y, nc, nm);
}

double get_shock_wave()
{
	double fvel = 446.11;
	double ftmp = temp[0];
	double vmp = sqrt(2.0*K*ftmp/mass[0]);
	double sc = fvel/vmp;
	double fs1 = sc + sqrt(sc*sc+2.0);
	double fs2 = 0.5*(1.0+sc*(2.0*sc-fs1));
	double qa = 3.0;
	if (sc < -3.0) {qa = fabs(sc)+1.0;}
	double u = 0, un = 0, a = 0;
	do{
		u = -qa+2.0*qa*frand();
		un = u + sc;
		if (u<0) {a = -1;continue;}
		a = (2.0*un/fs1)*exp(fs2-u*u);
 	} while (a < frand());
	return (un*vmp);
}

void rvelc(double *u, double *v, double vmp)
{
	double a = sqrt(-log(frand()));
	double b = 2.0*PI*frand();
	*u = a*sin(b)*vmp;
	*v = a*cos(b)*vmp;
}

void x_move(int krp, int x, int y, int nc, int nm)
{
	int zx;
//  	double fvel = 446.11;
// 	double dvel = 282.23;
	double vx;
	VECTOR newVX;
	double part1, part2, wt;
	PARTICLE *p, *b;
	p = &(camera[krp][x][y][nc][nm]);
	zx = ((p->v.x) > 0) ? (1) : (-1);
	copy_vector(&(p->v), &newVX);
	vx = fabs(p->v.x);
	wt = p->weight;
	part2 = (vx * dtp/SIZEX)*wt;
	part1 = (wt - part2);
	p->weight = part1;

	if (((KX - 1) == x) && (zx > 0)) 
	{
		zx = 0;
		newVX.x *= (-1.0);
		rvelc(&(newVX.y),&(newVX.z), sqrt(2.0*K*temp[0]/mass[0]));
	}
	if ((0 == x) && (zx < 0))
	{
		zx = 0;
 		newVX.x = get_shock_wave();
		rvelc(&(newVX.y),&(newVX.z), sqrt(2.0*K*temp[0]/mass[0]));
	}

	b = &(buf[(x+zx)][y][nm]);
 	b->weight += part2;
	if (frand() < (part2 / (b->weight)))
	{
		copy_vector(&newVX, &(b->v));
		b->mass = p->mass;
		b->diam = p->diam;
		b->rot_energy = p->rot_energy;
	}
}

void move_particles(int nc, int krp)
{
 	int x, nm;
	int y;
	for (x = 0; x < KX; x++)
	{
		for (y = 0; y < KY; y++)
		{
			for (nm = 0; nm < NM; nm++)
			{
				make_2d_move(krp, x, y, nc, nm);
			}
		}
	}
}

void move_camera(int nc, int krp)
{
	double curr_time = 0;
	while (curr_time < dtp) 
	{
		null_buf(nc, krp);
		move_particles(nc, krp);
		final_buf(nc, krp);
		curr_time += dtp_var[krp];
	}
}

