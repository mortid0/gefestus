#ifndef PICTURES_H
#define PICTURES_H

/* ������� ����� ���������� ������ ����������. ���� �� ����� ������� �������� - ���������� �������� 0 */

int SHOW_TEMP_FREQ;
int SHOW_CONC_FREQ;
int SHOW_MASS_FREQ;
#define SHOW_NALIP_FREQ 0
int SHOW_IMPULSE_FREQ;

/*! \fn void show_field(int num, double (*cell_func)(int, int, int), char *func_dir, int startx, int endx, int starty, int endy);
 * ������� �������� ���� ������������� �������� ������� ����(����������, ������������ � .�.)
 * \param num - ����� ���������� ����
 * \param cell_func(nc,x,y) - ������� ���� � �������� ����� ������(x,y) ��� ������������� ����������(nc)
 * \param func_dir - �������� ��������, � ������� ����� ����������� ������.
 * \param startx - ����� ������ �� �����������, � ������� ���������� ���� ���������
 * \param endx - ����� ������ �� �����������, �� ������� ������������� ���� ���������
 * \param starty - ����� ������ �� ���������, � ������� ���������� ���� ���������
 * \param endy - ����� ������ �� ���������, �� ������� ������������� ���� ���������
 */
extern void show_field(int num, double (*cell_func)(int, int, int), char *func_dir, int startx, int endx, int starty, int endy);

/*! ������� ���� ���������� ��� 2 �������. 1 ������� - ������� ����, 2 - �����.
 */
extern void show_temp_2flows(int num);

/*!������� ���� ������������ ��� 2 �������. 1 ������� - ������� ����, 2 - �����.
 */
extern void show_conc_2flows(int num);

/*!������� ���� vy ��� 2 �������.
 */
extern void show_vy_2flows(int num);

extern void show_vector(int num);

//extern void show_mass(int num);

//extern void show_nalip(int num);

/*!
 * ������� ����������� ��������� � ��������� �� �� �����.
 * \param step - ����� �������� ���������� ����.
*/
extern void print_pictures(int step);

#endif
