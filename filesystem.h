#ifndef FILESYSTEM_H
#define FILESYSTEM_H


/*!
 * ������� ��������� ��������� ��� ������ ���������.
 */
extern void make_dirs();

/*!
 * ��������� ������� ��������� ������ � ���������� �� ��������� ���� \param step � ���� �� ����.
 */
extern void save_camera(int step);

/*!
 * ��������� ��������� ������ � ���������� � ���������� ���� \param step
 */
extern void load_camera(int step);
#endif
