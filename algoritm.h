#ifndef ALGORITM_H
#define ALGORITM_H

/*! \fn void normirovka()
*	\brief ��������� �������� ������ � ������.
*
*	������� �������� �������� ������ ����� �������, ����� ����� ����� �� �������������� ���� ������ 1. ��� ��������� ������ ���� dtp ������ ����� �������� ������� ����������.
*/
void normirovka();

/*! \fn void mixture()
*	\brief ���������� ������������� ������ ����� ������������.
*
*	��������� ����� ������ ����� ������������, ��� ���� ������� �� ����� ����������
*	����� ������� � ������ ����������. ����� ��������� �������� ����� ���������� �����
*	���������� ������ ���� �� ����� ���� � ������ � ���������� � �������� ������������.
*/
void mixture();

/*!
*	\fn void inject(double time)
*	\brief ���������� ������������ ������ � ������ � ����� ������.
*
*	����������� � ������ ���������� ������� � ����� CONC_INJ � ������������ TEMP_INJ �
*	���� ������ CELL_INJ.
*
*/
void inject(double time);

void inject_left_right_half();

/*!
 *	\fn void inject_left(double time)
 *	\brief ���������� ������������ ������ � ����� ������ � ������.
 *
 *	����������� � ������ ���������� ������� � ����� CONC_INJ � ������������ TEMP_INJ �
 *	������ ������  ����� ������.
 *
 */
void inject_left(double time);	

/*!
 *	\fn void inject_right(double time)
 *	\brief ���������� ������������ ������ � ����� ������ � ������.
 *
 *	����������� � ������ ���������� ������� � ����� CONC_INJ � ������������ TEMP_INJ �
 *	������ ������  ����� ������.
 *
 */
void inject_right(double time);	

void change_plate(int krp); /* Smena plastin v slu4ae proleta */

void force(int krp);

void accurate_force(int krp);

#endif
