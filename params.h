#ifndef PARAMS_H
#define PARAMS_H

int KRP, NM, NC, KX, KY;

#define K 0.5 /* Boltzman constant*/

#define kDisp 1.0		/*  */
#define G_SLIP 0.0
#define SF 0.0			/* ������ ������� Qnew = (1-SF)*Qold */
#define EPS 0.0001

#define PI 3.1415926
double dtp;

//note that experiment time must be float, else look in main.c for time_steps
double EXPERIMENT_TIME;

double SIZEX;
double SIZEY;
double SIZEZ;
// double test_move1, test_move2;
// double nm_test_move1, nm_test_move2;

/*!
 * \struct VECTOR ��������� ������ � 3-������ ������������ � ������������ (x, y, z)
 */
typedef struct vector {
    double x, y, z;
} VECTOR;


/*!
 * \struct PARTICLE ��������� ������� � ��������� �����������
 * \param v - ��������
 * \param mass - �����
 * \param diam - �������
 * \param weight - ���(������������)
 */
typedef struct particle {
	double mass, diam, weight;
	VECTOR v;
	double rot_energy; // rotational internal energy
	int rot_degr; // number of rotational degrees of freedom
} PARTICLE;

/*! camera to work */
PARTICLE *****camera;//[KRP][KX][KY][NC][NM];

/*! plate array for storing some temp data, e.g. conc fields. */
double **temp_plate;//[KX][KY];

/*! storing component params */
double *conc, *temp, *diam, *mass;
int *num_rot_degrees, *num_vib_degrees;

/*! collision number */
int colcount;

/*! buf plate for moving */
PARTICLE ***buf;//[KX][KY][NM];

/*! buf plate for collision reaktions */
PARTICLE **buf_cell;//[NC][NM];

/*! array for sticking particles to walls. Used in moving.*/
double *kNalip;//[NC];

/* Friction in walls. Used in moving. */
#define kFric 1.0

//dirs for storing files with data
#define TEMP_DIR "temp2/temp"
#define CONC_DIR "temp2/conc"
#define IMPULSE_DIR "temp2/impulse"
#define data_file "temp2/data.txt"
#define MASS_DIR "temp2/mass"

#define SAVE_DIR "temp2/save"	

#define SAVE_FREQ 100		
#define ETALON_MASS (5.0E-26)		
#define PIC_MASS_COUNT 5
#define PIC_MASS_STEP 20	


#define NALIP_DIR "temp2/nalip"
#define NALIP_ETALON_MASS 1.0	/* etalonnaya massa dlya pods4eta nalipwih 4astic */
#define NALIP_MASS_COUNT 10	/* maksimalnaya massa 4astic, u4ityvaemyh pri nalipanii na steny */

double ***left_wall_nalip, ***right_wall_nalip;
double ***up_wall_nalip, ***down_wall_nalip; //[KX][NC][NALIP_MASS_COUNT];

#define MASTER_RANK 0		/* Index(rank) of master process */

double *dtp_var;

#define AMPUTATION_SIGNAL 16	

#define kProlet 0.5	
int k_max_prolet;	
int is_attention;
double conc_leave;	

double after_move, after_collision;

#include <stdio.h>

FILE *plot_file;
#define PLOT_FILE "temp2/plot.dat"
double parts[4];

#define SINGLE_MODE 1
#define MULTI_MODE 0
#define RUN_MODE SINGLE_MODE	
int my_rank; 
int num_proc;

/* ���� ������������ */
int IS_INJECT;		

 /* ������ ������������, ������������� �� ����� �� 0 �� TIME_INJ */
double CONC_INJ;
/* ����������� ����� */
double TEMP_INJ;	
/* ����� ����� */
double MASS_INJ;	
/* ������� ����� */
double DIAM_INJ;	

/* ����� ������, � ������� ������������� ����� */
int CELL_INJ;		
/* ����� ��������������  */
double TIME_INJ;	

/* ����������� ������ ������� ��� ������������. 0 - ��� ������, 1 - ������ ������ ������� */
//double kCollision;	

/*! ��������� ��� �������� ���������� �����.
	\param temp - ����������� �����
	\param kAcc - ����������� ����������� �����. ����������, ����� ����� ������ ���������� ��������, � ����� ���������. 0 - ���������� ���������, 1 - ���������.
	\param kFric - ����������� ������ � �����.
*/
typedef struct wall_params {
	double	temp, kAcc, kFriction;
} WALL_PARAMS;

WALL_PARAMS left_wall, right_wall, up_wall, down_wall;

//viscos - ����������� �������� ������� (w � �����)
#define viscos 0.75
double *svib;
double cell_temp, cell_density;
double *ctmp, *cdnn;
double alfa;
double elev[2][48];//level numbers start from 1 !!!
double sigma_three, colzv, avzv;
#define REF_VIB_TEMP 3395.0
#define Zrot 5.0
#define GAM 1.0
//6.67E-11
#define MAX_PLOT_VEL 4002
VECTOR **force_field;

double ***temp_field, ***conc_field, ***mass_field, **distr_x, **distr_y, ***impulse_field;

double *avr_vx, *avr_vy;
double step_dissoc, total_dissoc;
double fin_c[3], fin_m[3];

#endif
