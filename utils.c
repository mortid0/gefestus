#include <math.h>
#include "params.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <regex.h>
#include <string.h>
#include "utils.h"

// double frand()
// {
// 	int i = rand();
//  	return ran1(&i);
// // 	return (((double)rand())/(RAND_MAX+1.0));
// }

double mag2(struct vector *v)
{
	return ((v->x)*(v->x)+(v->y)*(v->y)+(v->z)*(v->z));
}

double mag(struct vector *v)
{
	return sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
}

void minus(struct vector *result, struct vector *v1, struct vector *v2)
{
	result->x = v1->x - v2->x;
	result->y = v1->y - v2->y;
	result->z = v1->z - v2->z;	
}

void plus(struct vector *result, struct vector *v1, struct vector *v2)
{
	result->x = v1->x + v2->x;
	result->y = v1->y + v2->y;
	result->z = v1->z + v2->z;
}

void multiply(struct vector *result, struct vector *v, double a)
{
	result->x = v->x * a;
	result->y = v->y * a;
	result->z = v->z * a;
}

void copy_vector(struct vector *from, struct vector *to)
{
	to->x = from->x;
	to->y = from->y;
	to->z = from->z;
}

double get_energy(PARTICLE *p)
{
	return (0.5*p->weight*p->mass*mag2(&p->v));
}

void copy_particle(struct particle *from, struct particle *to)
{
	to->rot_energy = from->rot_energy;
	to->rot_degr = from->rot_degr;
	to->weight = from->weight;
	to->mass = from->mass;
	to->diam = from->diam;
	to->v.x = from->v.x;
	to->v.y = from->v.y;
	to->v.z = from->v.z;
}

double get_temp(PARTICLE *p, VECTOR *v)
{
	return p->weight*p->mass*(pow(p->v.x - v->x, 2)+pow(p->v.y - v->y, 2)+pow(p->v.z - v->z, 2))/(K*3);
}

double get_t(int nc, int x, int y)
{
	int krp, nm;
	double result = 0, norm = 0;
	VECTOR v0;
// 	PARTICLE *p;
	v0.x = v0.y = v0.z = 0;
// 	for (krp = 0; krp < KRP; krp++)
// 	{
// 		for (nm = 0; nm < NM; nm++)
// 		{
// 			p = &camera[krp][x][y][nc][nm];
// 			v0.x += p->v.x * p->weight;
// 			v0.y += p->v.y * p->weight;
// 			v0.z += p->v.z * p->weight;
// 			norm += p->weight;
// 		}
// 	}
// 	multiply(&v0, &v0, (1.0 / norm));
      
	for (krp = 0; krp < KRP; krp++)
	{
		for (nm = 0; nm < NM; nm++)
		{
			result += get_energy(&camera[krp][x][y][nc][nm]);
			norm += camera[krp][x][y][nc][nm].weight;
		}
	}
	if (0 == norm) return 0;
	return (result/norm);
// 	return (((result/norm)-293.0)/(407.81 - 293.0));
}

double get_xt(int nc, int x, int y)
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
			p = &camera[krp][x][y][nc][nm];
			result += p->weight*p->mass*(p->v.x - v0.x)*(p->v.x - v0.x)/(3*K);
		}
	}
	return (result/norm);
}

double get_yt(int nc, int x, int y)
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
			p = &camera[krp][x][y][nc][nm];
			result += p->weight*p->mass*(p->v.x - v0.x)*(p->v.x - v0.x)/(3*K);
		}
	}
	return (result/norm);
}

double get_qflow(int nc, int x, int y)
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
			p = &camera[krp][x][y][nc][nm];
			result += get_temp(p, &v0) * (p->v.x - v0.x) * (1.5 * K);
		}
	}
	return (result/norm);
}


int match(char match[N_MATCHES][MATCH_LEN], char *pattern, char *string)
{
	regex_t regex;
	regmatch_t pmatch[N_MATCHES];
	size_t re, i;
	if (0 != (re = regcomp(&regex, pattern, REG_EXTENDED)))
	{
		printf("invalid expression: %s\n",pattern);
		return 1;
	}
	re = regexec(&regex, string, N_MATCHES, pmatch, 0);
	for (i=0; !re && i<=regex.re_nsub ; i++)
	{
		strncpy(match[i], string + pmatch[i].rm_so, pmatch[i].rm_eo - pmatch[i].rm_so);
		match[i][pmatch[i].rm_eo-pmatch[i].rm_so]='\0';
	}
	regfree(&regex);
	return 0;
}

double get_hi2(double m, double d)
{
	double a,v;
	do
	{
		v = -3.0 + 6.0 * frand();
		a = exp(-v*v);
	} while ( a <frand());
	return (v*d+m);
}

void get_vy_half_maxwell(VECTOR *v, double mass, double temp)
{
	double d = sqrt(temp/mass);
	v->y = d*sqrt(-log(frand()));
	double vr = d*sqrt(-log(frand()));
	double phi = 2.0*PI*frand();
	v->z = vr*sin(phi);
	v->x = vr*cos(phi);
}

void get_vx_half_maxwell(VECTOR *v, double mass, double temp)
{
	double d = sqrt(temp/mass);
	v->x = d*sqrt(-log(frand()));
	double vr = d*sqrt(-log(frand()));
	double phi = 2.0*PI*frand();
	v->y = vr*sin(phi);
	v->z = vr*cos(phi);
}

void move_particle_to_buf(PARTICLE *p, PARTICLE *b, VECTOR *v, double new_weight)
{
	b->weight += new_weight;
	if (frand() > (new_weight / (b->weight))) {return;}

	b->v.x = v->x;
	b->v.y = v->y;
	b->v.z = v->z;
	b->diam = p->diam;
	b->mass = p->mass;
	b->rot_energy = p->rot_energy;
	b->rot_degr = p->rot_degr;
}

double gam(double x)
{
	double a = 1.0, y = x;
	if (y < 1.0) {a = a/y;}
	else {
	    y -= 1;
	    while (y>1){
		a *= y;
		y -= 1;
	    }
	}
	return (a*(1.0-0.5748646*y+0.9512363*y*y - 0.6998688*y*y*y + 0.42455499*y*y*y*y - 0.1010678*y*y*y*y*y));
}

//returns Larsen-Borgnakke internal energy redistribution
double lbs(double XMA, double XMB)
{
    double erm, p;
    do {
	erm = frand();
	if ((XMA < 0.000001) || (XMB < 0.000001))
	{
	    if ((XMA < 0.000001) && (XMB < 0.000001)) {return erm;}
	    if (XMA < 0.000001) {p = pow(1.0 - erm, XMB);}
	    if (XMB < 0.000001) {p = pow(1.0 - erm, XMA);}
	} else {p = pow((XMA+XMB)*erm/XMA,XMA)*pow((XMA+XMB)*(1.0-erm)/XMB,XMB);}
    }while (frand() > p);
    return erm;
}

