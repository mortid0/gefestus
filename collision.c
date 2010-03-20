#include "collision.h"
#include "params.h"
#include <math.h>
#include "utils.h"

#define Qdiss (1.0E20)
#define ecd 53.0*K
#define kFricCol 1.0

#define TFAC 0.1
#define DFI 2.0
#define sigma_grav 1.0


void elastic_collision(int r)
{
	int ix, iy, l, i, j, lm;
	int lp, mp;
	struct vector vrc, vccm;
	PARTICLE *pl, *pm;
	double sigma, a, b, g;
	double sm, rmm, rml;
	for (iy = 0; iy < KY; iy++)
	{
		for (ix = 0; ix < KX; ix++)
		{
			for (l = 0; l < NC; l++)
			{
				for (i = 0; i < NM - 1; i++)
				{
					for (j = i + 1; j < NM; j++)
					{
						pm = &camera[r][ix][iy][l][i];
						pl = &camera[r][ix][iy][l][j];
// 						sigma = PI * pow(pl->diam+pm->diam, 2)/4.0;
						g = sqrt(pow(pl->v.x - pm->v.x,2) + pow(pl->v.y - pm->v.y,2) + pow(pl->v.z - pm->v.z,2));
						double rm = pl->mass*pm->mass/(pl->mass+pm->mass);
						sigma = PI * 1.5 * 1.5 * pow(1.0/(rm*g*g),0.25);
						if (pl->weight > pm->weight)
						{
							a = pl->weight * sigma * g * dtp_var[r] * NM / (NM - 1.0);
						}
						else
						{
							a = pm->weight * sigma * g * dtp_var[r] * NM / (NM - 1.0);
						}
						b = frand();
						if (b < a)
						{
							double wa = sqrt(frand()) * 1.5;
							double eps = 2.0* PI * frand();
							a = wa*(1.26233+wa*(1.84145+wa*(-87881+wa*(20.3313 + wa*(-23.8155+wa*(14.5046+wa*(-4.42027+wa*0.535193)))))));
							double chi = PI - 2.0*a;
							double cc = cos(chi);	double sc = sin(chi);
							double ce = cos(eps);	double se = sin(eps);
							VECTOR vrc_old;
							vrc_old.x = pl->v.x - pm->v.x;vrc_old.y = pl->v.y - pm->v.y;vrc_old.z = pl->v.z - pm->v.z;
							
							vrc.x = cc * vrc_old.x + sc*se*sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
							vrc.y = cc * vrc_old.y + sc*(g*vrc_old.z*ce-vrc_old.x*vrc_old.y*se)/sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
							vrc.z = cc * vrc_old.y - sc*(g*vrc_old.y*ce+vrc_old.x*vrc_old.z*se)/sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
							
							lp = 1;
							mp = 1;
							b = frand();
							if (pm->weight > EPS)
							{
								a = pl->weight / pm->weight;
							}
							else
							{
								a = 1.0/EPS;
							}
							if (b > a) {mp = 0;}
							else if (b > (1.0/a)) {lp = 0;}
							vccm.x = (pl->v.x + pm->v.x) / 2.0;
							vccm.y = (pl->v.y + pm->v.y) / 2.0;
							vccm.z = (pl->v.z + pm->v.z) / 2.0;
							if (1 == lp)
							{
								pl->v.x = (vccm.x + vrc.x / 2.0);
								pl->v.y = (vccm.y + vrc.y / 2.0);
								pl->v.z = (vccm.z + vrc.z / 2.0);
							}
							if (1 == mp)
							{
								pm->v.x = (vccm.x - vrc.x / 2.0);
								pm->v.y = (vccm.y - vrc.y / 2.0);
								pm->v.z = (vccm.z - vrc.z / 2.0);
							}
						}
					}
				}
				if (l < NC-1)
				{
					for (lm = l; lm < NC; lm++)
					{
			
						for (i = 0; i < NM; i++)
						{
							pl = &camera[r][ix][iy][l][i];
							for (j = 0; j < NM; j++)
							{
								pm = &camera[r][ix][iy][lm][j];
								
								g = sqrt(pow(pl->v.x - pm->v.x,2) + pow(pl->v.y - pm->v.y,2) + pow(pl->v.z - pm->v.z,2));
// 								sigma = PI * pow(pl->diam + pm->diam, 2) / 4.0;
								double rm = pl->mass*pm->mass/(pl->mass+pm->mass);
								sigma = PI * 1.5 * 1.5 * pow(1.0/(rm*g*g),0.25);
								if (pl->weight > pm->weight)
								{
									a = pl->weight * sigma * g * dtp_var[r];
								}
								else 
								{
									a = pm->weight * sigma * g * dtp_var[r];
								}
								b = frand();
								if (b < a)
								{
									double wa = sqrt(frand()) * 1.5;
									double eps = 2.0* PI * frand();
									a = wa*(1.26233+wa*(1.84145+wa*(-87881+wa*(20.3313 + wa*(-23.8155+wa*(14.5046+wa*(-4.42027+wa*0.535193)))))));
									double chi = PI - 2.0*a;
									double cc = cos(chi);	double sc = sin(chi);
									double ce = cos(eps);	double se = sin(eps);
									VECTOR vrc_old;
									vrc_old.x = pl->v.x - pm->v.x;vrc_old.y = pl->v.y - pm->v.y;vrc_old.z = pl->v.z - pm->v.z;
									
									vrc.x = cc * vrc_old.x + sc*se*sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
									vrc.y = cc * vrc_old.y + sc*(g*vrc_old.z*ce-vrc_old.x*vrc_old.y*se)/sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
									vrc.z = cc * vrc_old.y - sc*(g*vrc_old.y*ce+vrc_old.x*vrc_old.z*se)/sqrt(vrc_old.y*vrc_old.y+vrc_old.z*vrc_old.z);
									
									lp = 1;
									mp = 1;
									sm = pl->mass + pm->mass;
									rml = pl->mass / sm;
									rmm = pm->mass / sm;
									b = frand();
									if (pm->weight > EPS)
									{
										a = pl->weight / pm->weight;
									}
									else
									{
										a = 1.0 / EPS;
									}
									if (b > a) {mp = 0;}
									else if (b > (1.0/a)) {lp = 0;}
												
									vccm.x = rml * pl->v.x + rmm * pm->v.x;
									vccm.y = rml * pl->v.y + rmm * pm->v.y;
									vccm.z = rml * pl->v.z + rmm * pm->v.z;
									
									if (1 == lp)
									{
										pl->v.x = vccm.x + vrc.x * rmm;
										pl->v.y = vccm.y + vrc.y * rmm;
										pl->v.z = vccm.z + vrc.z * rmm;
									}
									if (1 == mp)
									{
										pm->v.x = vccm.x - vrc.x * rml;
										pm->v.y = vccm.y - vrc.y * rml;
										pm->v.z = vccm.z - vrc.z * rml;
									}	
								}
							}
						}
					}
				}
			}
		}
	}
}


double collision_probability(PARTICLE *pl, PARTICLE *pm)
{
	if ((0==pl->weight)||(0==pm->weight)) {return 0;}
	double sigma, g, a = (pl->weight > pm->weight)?pl->weight:pm->weight;
	g = sqrt(pow(pl->v.x - pm->v.x, 2) + pow(pl->v.z - pm->v.z, 2) + pow(pl->v.z - pm->v.z, 2));

	sigma = 0.25 * PI * pow(pl->diam+pm->diam, 2) * sigma_grav * sigma_grav;
// 	double cref = 2.0 * sqrt(2.0*K*273.0/(rmass*PI));
// 	sigma *= pow(cref/g, 2.0 * (viscos - 0.5));
	a *= sigma * g * dtp;
	
	if (a>1) {printf("a = %e;w = %e; g = %e; sigma = %e; pl.diam = %e; pm.diam = %e\n", a, (pl->weight>pm->weight)?(pl->weight):(pm->weight), g, sigma, pl->diam, pm->diam);}
	return a;
}

//return part of collision that can redistribute rotational energy
double get_zrot2(PARTICLE *tmp_particle)
{
// 	if (tmp_particle->vib_level > MLEV) {return 0;}
	if (2 == tmp_particle->rot_degr) {return 0.2;}
	if (3 == tmp_particle->rot_degr) {return 0.1;}
	return 0;
}

void redistribute_rot_energy2(PARTICLE *tmp_particle, double *Et)
{
	double Ecc = (*Et) + tmp_particle->rot_energy;
	double Eci = tmp_particle->rot_energy;
	double Erm, XIB = 2.5 - viscos;
	if (2 == tmp_particle->rot_degr){Erm = 1.0 - pow(frand(), (1.0/XIB));}
	else {double XIA = tmp_particle->rot_degr*0.5;Erm = lbs((XIA - 1.0), (XIB - 1.0) );}
	tmp_particle->rot_energy = Erm * Ecc;
	(*Et) += (Eci - tmp_particle->rot_energy);
}

void change_internal_energy_after_collision2(PARTICLE *pl, PARTICLE *pm, VECTOR *vrc, double *vr, double *vrr)
{
	double rmass = (pl->mass/(pl->mass+pm->mass)) * pm->mass;
	double Et = 0.5 * rmass * (*vrr);
	int irt = 0;
	int i;
	PARTICLE *tmp_particle;
	for (i = 0; i < 2; i++)
	{
		tmp_particle = pl;
		if (1 == i) {tmp_particle = pm;}
		if ((tmp_particle->rot_degr > 0) && (frand() < get_zrot2(tmp_particle)))
		{
			irt = 1;
			redistribute_rot_energy2(tmp_particle, &Et);
		}
	}
	if (1 == irt)
	{
		double a = sqrt(2.0 * Et / rmass);
		vrc->x *= (a/(*vr));
		vrc->y *= (a/(*vr));
		vrc->z *= (a/(*vr));
		(*vr) = a;
	}
}

void change_velocity_after_collision2(PARTICLE *pl, PARTICLE *pm, VECTOR *vrc, double *vr)
{
	VECTOR vccm;
	double a, b;
	double sm, rmm, rml, rmass;

	sm = pl->mass + pm->mass;
	rml = pl->mass / sm;
	rmm = pm->mass / sm;
	rmass = rml * pm->mass;
	if ((*vr) > sqrt(2.0*Qdiss/rmass)) {(*vr) = sqrt(pow((*vr),2) - 2.0*Qdiss/rmass);}
	else {(*vr) = 0;}
	vccm.x = rml * pl->v.x + rmm * pm->v.x;
	vccm.y = rml * pl->v.y + rmm * pm->v.y;
	vccm.z = rml * pl->v.z + rmm * pm->v.z;
	
// 	double spm4 = 0.605;
	b = 2.0 * frand() - 1.0;
	a = sqrt(1.0 - b*b);
	VECTOR vrcp;
	vrcp.x = b * (*vr);
	double c = 2.0 * PI * frand();
	vrcp.y = a * cos(c) * (*vr);
	vrcp.z = a * sin(c) * (*vr);
	
	pl->v.x = vccm.x + vrcp.x * rmm;
	pl->v.y = vccm.y + vrcp.y * rmm;
	pl->v.z = vccm.z + vrcp.z * rmm;

	pm->v.x = vccm.x - vrcp.x * rml;
	pm->v.y = vccm.y - vrcp.y * rml;
	pm->v.z = vccm.z - vrcp.z * rml;
}

void change_velocity_after_gravitation(PARTICLE *pl, PARTICLE *pm, VECTOR *vrc, double *vr, double b)
{
	VECTOR vccm;
	double chi, psi, wid;
	double sm, rmm, rml, rmass;

	sm = pl->mass + pm->mass;
	rml = pl->mass / sm;
	rmm = pm->mass / sm;
	rmass = rml * pm->mass;
	wid = b * 0.5 * (pl->diam + pm->diam);

	vccm.x = rml * pl->v.x + rmm * pm->v.x;
	vccm.y = rml * pl->v.y + rmm * pm->v.y;
	vccm.z = rml * pl->v.z + rmm * pm->v.z;
	
	chi = -2.0 * atan(GAM * sm / ((*vr)*(*vr) * wid));
	
	psi = 2.0 * PI * frand();
      
	VECTOR vrcp;
	vrcp.x = cos(chi) * vrc->x + sin(chi)*sin(psi)*sqrt(vrc->y*vrc->y + vrc->z*vrc->z);
	vrcp.y = cos(chi) * vrc->y + sin(chi)*((*vr)*vrc->z*cos(psi) - vrc->x*vrc->y*sin(psi))/sqrt(vrc->y*vrc->y + vrc->z*vrc->z);
	vrcp.z = cos(chi) * vrc->z - sin(chi)*((*vr)*vrc->y*cos(psi) + vrc->x*vrc->z*sin(psi))/sqrt(vrc->y*vrc->y + vrc->z*vrc->z);
	
	pl->v.x = vccm.x + vrcp.x * rmm;
	pl->v.y = vccm.y + vrcp.y * rmm;
	pl->v.z = vccm.z + vrcp.z * rmm;

	pm->v.x = vccm.x - vrcp.x * rml;
	pm->v.y = vccm.y - vrcp.y * rml;
	pm->v.z = vccm.z - vrcp.z * rml;
}


void weight_method2(PARTICLE *pl, PARTICLE *new_pl, PARTICLE *pm, PARTICLE *new_pm)
{
	int lp = 1, mp = 1;
	double a = (new_pm->weight > 0)?(new_pl->weight / new_pm->weight):(1E6);
	
	double b = frand();
	if (b > a) {mp = 0;}
	else if (b > (1.0/a)) {lp = 0;}
	if (1 == lp) {copy_particle(new_pl, pl);}
	if (1 == mp) {copy_particle(new_pm, pm);}
}

double make_cluster(PARTICLE *new_pl, PARTICLE *new_pm, double vr)
{
	double rmass = (new_pl->mass/(new_pl->mass+new_pm->mass))*new_pm->mass;
 	double cref = temp[0];
	double clust_prob = exp(-10.0*vr/cref);
	if (clust_prob < 0.005) return 1.0;
	double cluster_weight = clust_prob * ((new_pl->weight>new_pm->weight)?(new_pm->weight):(new_pl->weight));
	new_pl->weight -= cluster_weight;
	new_pm->weight -= cluster_weight;
	PARTICLE cluster;
	cluster.weight = cluster_weight;
	cluster.mass = new_pl->mass + new_pm->mass;
	cluster.diam = pow(pow(new_pl->diam,3) + pow(new_pm->diam,3),0.33);
	VECTOR cluster_vel;
	cluster_vel.x = (new_pl->mass*new_pl->v.x + new_pm->mass*new_pm->v.x)/(new_pl->mass+new_pm->mass);
	cluster_vel.y = (new_pl->mass*new_pl->v.y + new_pm->mass*new_pm->v.y)/(new_pl->mass+new_pm->mass);
	cluster_vel.z = (new_pl->mass*new_pl->v.z + new_pm->mass*new_pm->v.z)/(new_pl->mass+new_pm->mass);
	double mag_vel = mag(&cluster_vel);
	double part_vel_to_remain = 1.0;
	if (mag_vel > sqrt(2.0*Qdiss/rmass)) {part_vel_to_remain = sqrt(pow(mag_vel,2) - 2.0*Qdiss/rmass)/mag_vel;}
	else {part_vel_to_remain = 0.0;}
 	multiply(&cluster_vel, &cluster_vel, part_vel_to_remain);
	copy_vector(&cluster_vel, &(cluster.v));
	int clust_nm = rand() % NM;
	PARTICLE *buff = &(buf_cell[2][clust_nm]);
	move_particle_to_buf(&cluster, buff, &(cluster.v), cluster.weight);
	return part_vel_to_remain;
// 	printf("Cluster mass=%e; weight = %e\n", cluster.mass, cluster.weight);
}

void manage_collision(PARTICLE *pl, PARTICLE *pm)
{
	VECTOR vrc;

	minus(&vrc, &(pl->v), &(pm->v));
	double vr = mag(&vrc);
	
// 	double vrr = vr*vr;
	PARTICLE new_pl, new_pm;
	copy_particle(pl, &new_pl);
	copy_particle(pm, &new_pm);
//	double b = sqrt(frand())*sigma_grav;
//	if (b <= 1.0) // b < (pl->diam + pm->diam) /2 - happens collision
//	{
// 	change_internal_energy_after_collision2(&new_pl, &new_pm, &vrc, &vr, &vrr);
	  vr *= make_cluster(&new_pl, &new_pm, vr);
	  change_velocity_after_collision2(&new_pl, &new_pm, &vrc, &vr);
//	}
//	else // happens force interaction
//	{
//	  change_velocity_after_gravitation(&new_pl, &new_pm, &vrc, &vr, b);
//	}
	  pl->weight = new_pl.weight;
	  pm->weight = new_pm.weight;
	  weight_method2(pl, &new_pl, pm, &new_pm);
}

void internal_component_collision(PARTICLE *pl, PARTICLE *pm)
{
// 	printf("collision prob = %e\n",collision_probability(pl, pm));
	if (frand() > collision_probability(pl, pm) * NM / (NM - 1.0)) return;
	
	manage_collision(pl, pm);
}

void external_component_collision(PARTICLE *pl, PARTICLE *pm)
{
	if (frand() > collision_probability(pl, pm)) return;

	manage_collision(pl, pm);
}


void inelastic_collision(int r)
{
	int ix, iy, l, i, j, lm;

	PARTICLE *pa, *pb;
//	printf("col_prob = %e\n", collision_probability(&camera[0][0][0][0][0], &camera[0][0][0][0][1]));
	for (iy = 0; iy < KY; iy++)
	{
		for (ix = 0; ix < KX; ix++)
		{
			for (l = 0; l < NC; l++)
			{
				for (i = 0; i < NM; i++)
				{
					buf_cell[l][i].weight = 0;
				}
			}
			for (l = 0; l < NC; l++)
			{
				for (i = 0; i < NM - 1; i++)
				{
					for (j = i+1; j < NM; j++)
					{
						pa = &camera[r][ix][iy][l][i];
						pb = &camera[r][ix][iy][l][j];
						internal_component_collision(pa, pb);
					}
				}
				if (l > 0)
				{
					for (lm = 0; lm < l; lm++)
					{
						for (i = 0; i < NM; i++)
						{
							for (j = 0; j < NM; j++)
							{
								pa = &camera[r][ix][iy][l][i];
								pb = &camera[r][ix][iy][lm][j];
								external_component_collision(pa, pb);
							}
						}
					}
				}
			}
			for (l = 0; l < NC; l++)
			{
				for (i = 0; i < NM; i++)
				{
					pa = &camera[r][ix][iy][l][i];
					pb = &buf_cell[l][i];
					move_particle_to_buf(pb, pa, &(pb->v), pb->weight);
				}
			}
		}
	}
}
