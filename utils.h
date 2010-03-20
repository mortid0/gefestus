/*!	\file utils.h
	\brief ����, �������� �������� ��������������� �������

	� ������ ����� ���������� ������� ��� ������ �� ���������� �������, ��������� � ����������� �����������.
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "params.h"

/*!	\def frand()
	\brief ������, ������������ ��������� ��������.

	������ ���������� ��������� ����� ���� double ���������� �������������� �� 0 �� 1.
*/

#define frand() (((double)rand())/(RAND_MAX+1.0))

//extern double frand();
/*!	\def randomize()
	\brief �������� ����� ������������������ ��������� �����.

	������ ������������� ����� ��������� �������� ��� ��������� ��c��������������� �����.
*/

#define randomize() (srand((unsigned int)time((time_t *)NULL)))

/*!	\def N_MATCHES

	������ ������������� ���������� ������, ����������� � ������, � ���������� ���������.
*/

#define N_MATCHES 10

/*!	\def MATCH_LEN

	������������� ����� ������ ���������� ����������� ���������.
*/

#define MATCH_LEN 132

/*!	\fn double mag2(VECTOR *v)

	\brief ���������� ����� ������� �������� v � ��������.
	\param v ��������� �� ������ ��������.
*/

extern double mag2(VECTOR *v);

/*!	\fn double mag(VECTOR *v)

	\brief ���������� ����� ������� �������� v.
	\param v ��������� �� ������ ��������.
*/

extern double mag(VECTOR *v);

extern void minus(VECTOR *result, VECTOR *v1, VECTOR *v2);

extern void plus(VECTOR *result, VECTOR *v1, VECTOR *v2);

extern void multiply(VECTOR *result, VECTOR *v, double a);

extern double get_energy(PARTICLE *p);

/*!
 *	\fn double get_t(int nc, int x, int y);
 * 	\brief ���������� ����������� � ������
 * 	\param nc ����� ����������
 * 	\param x ���������� �� �����������
 * 	\param y ���������� �� ���������
 *
 * 	������� ����������� ��������� ���������� nc � ������ ������ (x, y).
 * 	� �������� ����������� ������� �������� �������������� ������� ���
 * 	����� ����� �������� ������ � ������ ������.
 * 
 */
extern double get_t(int nc, int x, int y);

/*!
 *	\fn double get_xt(int nc, int x, int y);
 * 	\brief ���������� ����������(x) ����������� � ������
 * 	\param nc ����� ����������
 * 	\param x ���������� �� �����������
 * 	\param y ���������� �� ���������
 *
 * 	������� ����������� ��������� ���������� nc �� ��� X � ������ ������ (x, y).
 * 	� �������� ����������� ������� �������� �������������� ������� ���
 * 	����� ����� �������� ������ � ������ ������.
 * 
 */
extern double get_xt(int nc, int x, int y);

/*!
 *	\fn double get_yt(int nc, int x, int y);
 * 	\brief ���������� ����������(y) ����������� � ������
 * 	\param nc ����� ����������
 * 	\param x ���������� �� �����������
 * 	\param y ���������� �� ���������
 *
 * 	������� ����������� ��������� ���������� nc �� ��� Y � ������ ������ (x, y).
 * 	� �������� ����������� ������� �������� �������������� ������� ���
 * 	����� ����� �������� ������ � ������ ������.
 * 
 */
extern double get_yt(int nc, int x, int y);

/*!
 *	\fn double get_qflow(int nc, int x, int y);
 * 	\brief ���������� ����� ����� � ������ �� ��� �
 * 	\param nc ����� ����������
 * 	\param x ���������� �� �����������
 * 	\param y ���������� �� ���������
 *
 * 	������� �������� ����� ��������� ���������� nc � ������ ������ (x, y).
 * 	� �������� ����� ������� ������������ �������������� ������� ���
 * 	����� ����� �������� ������ � ������ ������.
 *
 */
extern double get_qflow(int nc, int x, int y);


/*!
 *	\fn double get_temp(PARTICLE *p);
 * 	\brief ���������� ����������� �������
 * 	\param *p ��������� �� �������
 *
 * 	������� ����������� ������� *p. � ������� ����������� ������������
 * 	�������������� ������� �������.
 * 	
 */
extern double get_temp(PARTICLE *p, VECTOR *v);

extern void copy_particle(struct particle *from, struct particle *to);

extern void copy_vector(struct vector *from, struct vector *to);

/*! \fn int match(char match[N_MATCHES][MATCH_LEN], char *pattern, char *string)
 * \brief ����� ����������� ��������� � ������
 * \param pattern - ���������� ���������
 * \param string - ������ � ������� ���������� �����
 * \param match[][] - ������, � ������� ����������� ���������� ������. match[i] �������� ������� ������� � i-�� ������ ����������� ���������.
 * 
 * ������� ���� ������� � ������ � ���������� ��������� �������� � ���� ������� ��������.
 * 0 ������� - ��� ��������� ���������, ��������� - ��������������� ��������� � �������
 * ������� (��. ���������� ���������). ���� ���������� �� �������, �� ������������ NULL.
 * ��� ���������� ����� ���������� ���������� -lc -lcompat.
 * ����� ������������� ������������ �������� ������� ����������.
 */
extern int match(char match[N_MATCHES][MATCH_LEN], char *pattern, char *string);


/*! \fn double get_hi2(double m, double d)
 *	������� ���������� ��������� ��������, ������� hi^2 �������������. ��� �������� ����� m, � ��������� d. �������� ����� � ��������� �� -3d �� 3d.
 */
extern double get_hi2(double m, double d);

/*! \fn void get_vy_half_maxwell(VECTOR *v, double mass, double temp)
 * 	���������� ������ v � ������������ ��������������� � ���������������� (y > 0) � ������������ temp � ������ mass.
 */
extern void get_vy_half_maxwell(VECTOR *v, double mass, double temp);


/*! \fn void get_vx_half_maxwell(VECTOR *v, double mass, double temp)
 * 	���������� ������ v � ������������ ��������������� � ���������������� (x > 0) � ������������ temp � ������ mass.
 */
extern void get_vx_half_maxwell(VECTOR *v, double mass, double temp);

/*!
 * \fn void move_particle_to_buf(PARTICLE *p, PARTICLE *b, VECTOR *v, double new_weight)
 * 	����������� ��� ������� b �� new_weight � � ������������ new_weight / b->weight ������� � b ������� � ������, ��������� � ��������� p->mass, p->diam � v ��������������.
 */
extern void move_particle_to_buf(PARTICLE *p, PARTICLE *b, VECTOR *v, double new_weight);
/*!
 * \fn double gam(double x)
 * ���������� ����� ������� �� \param x ������������ 0
 */
extern double gam(double x);

/*!
 * \fn double lbs(double XIA, double XIB)
 * ���������� ������������� ���������� ������� 
 */
extern double lbs(double XIA, double XIB);

#endif
