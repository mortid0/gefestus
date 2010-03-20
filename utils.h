/*!	\file utils.h
	\brief Файл, содержит полезные вспомогательные функции

	В данном файле содержатся функции для работы со случайными числами, векторами и регулярными выражениями.
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "params.h"

/*!	\def frand()
	\brief Макрос, возвращающий случайное значение.

	Макрос возвращает случайное число типа double равномерно распределенное от 0 до 1.
*/

#define frand() (((double)rand())/(RAND_MAX+1.0))

//extern double frand();
/*!	\def randomize()
	\brief Начинает новую последовательность случайных чисел.

	Макрос устанавливает новое начальное значение для случайной поcледовательности чисел.
*/

#define randomize() (srand((unsigned int)time((time_t *)NULL)))

/*!	\def N_MATCHES

	Макрос устанавливает количество атомов, заключенных в скобки, в регулярном выражении.
*/

#define N_MATCHES 10

/*!	\def MATCH_LEN

	Устанавливает длину строки найденного регулярного выражения.
*/

#define MATCH_LEN 132

/*!	\fn double mag2(VECTOR *v)

	\brief Возвращает длину вектора скорости v в квадрате.
	\param v указатель на вектор скорости.
*/

extern double mag2(VECTOR *v);

/*!	\fn double mag(VECTOR *v)

	\brief Возвращает длину вектора скорости v.
	\param v указатель на вектор скорости.
*/

extern double mag(VECTOR *v);

extern void minus(VECTOR *result, VECTOR *v1, VECTOR *v2);

extern void plus(VECTOR *result, VECTOR *v1, VECTOR *v2);

extern void multiply(VECTOR *result, VECTOR *v, double a);

extern double get_energy(PARTICLE *p);

/*!
 *	\fn double get_t(int nc, int x, int y);
 * 	\brief Определяет температуру в ячейке
 * 	\param nc номер компонента
 * 	\param x координата по горизонтали
 * 	\param y координата по вертикали
 *
 * 	Снимает температуру заданного компонента nc в данной ячейке (x, y).
 * 	В качестве температуры берется тепловая поступательная энергия без
 * 	учета общей скорости потока в данной ячейке.
 * 
 */
extern double get_t(int nc, int x, int y);

/*!
 *	\fn double get_xt(int nc, int x, int y);
 * 	\brief Определяет продольную(x) температуру в ячейке
 * 	\param nc номер компонента
 * 	\param x координата по горизонтали
 * 	\param y координата по вертикали
 *
 * 	Снимает температуру заданного компонента nc по оси X в данной ячейке (x, y).
 * 	В качестве температуры берется тепловая поступательная энергия без
 * 	учета общей скорости потока в данной ячейке.
 * 
 */
extern double get_xt(int nc, int x, int y);

/*!
 *	\fn double get_yt(int nc, int x, int y);
 * 	\brief Определяет продольную(y) температуру в ячейке
 * 	\param nc номер компонента
 * 	\param x координата по горизонтали
 * 	\param y координата по вертикали
 *
 * 	Снимает температуру заданного компонента nc по оси Y в данной ячейке (x, y).
 * 	В качестве температуры берется тепловая поступательная энергия без
 * 	учета общей скорости потока в данной ячейке.
 * 
 */
extern double get_yt(int nc, int x, int y);

/*!
 *	\fn double get_qflow(int nc, int x, int y);
 * 	\brief Определяет поток тепла в ячейке по оси Х
 * 	\param nc номер компонента
 * 	\param x координата по горизонтали
 * 	\param y координата по вертикали
 *
 * 	Снимает тепловой поток заданного компонента nc в данной ячейке (x, y).
 * 	В качестве тепла берется кинетическая поступательная энергия без
 * 	учета общей скорости потока в данной ячейке.
 *
 */
extern double get_qflow(int nc, int x, int y);


/*!
 *	\fn double get_temp(PARTICLE *p);
 * 	\brief Определяет температуру частицы
 * 	\param *p указатель на частицу
 *
 * 	Снимает температуру частицы *p. В качеств температуры используется
 * 	поступательная энергия частицы.
 * 	
 */
extern double get_temp(PARTICLE *p, VECTOR *v);

extern void copy_particle(struct particle *from, struct particle *to);

extern void copy_vector(struct vector *from, struct vector *to);

/*! \fn int match(char match[N_MATCHES][MATCH_LEN], char *pattern, char *string)
 * \brief Поиск регулярного выражения в строке
 * \param pattern - регулярное выражение
 * \param string - строка в которой происходит поиск
 * \param match[][] - массив, в который сохраняются результаты поиска. match[i] содержит строчку стоящую в i-ой скобке регулярного выражения.
 * 
 * Функция ищет паттерн в строке и возвращает найденное значение в виде массива подстрок.
 * 0 элемент - это найденное выражение, следующие - соответствующие выражения в круглых
 * скобках (см. регулярные выражения). Если совпадение не найдено, то возвращается NULL.
 * Для компиляции нужно подключить библиотеки -lc -lcompat.
 * После использования возвращаемое значение следует освободить.
 */
extern int match(char match[N_MATCHES][MATCH_LEN], char *pattern, char *string);


/*! \fn double get_hi2(double m, double d)
 *	Функция возвращает случайное значение, имеющее hi^2 распределение. Мат ожидание равно m, а дисперсия d. Значения лежат в диапазоне от -3d до 3d.
 */
extern double get_hi2(double m, double d);

/*! \fn void get_vy_half_maxwell(VECTOR *v, double mass, double temp)
 * 	Возвращает вектор v с компонентами распределенными в полупространство (y > 0) с температурой temp и массой mass.
 */
extern void get_vy_half_maxwell(VECTOR *v, double mass, double temp);


/*! \fn void get_vx_half_maxwell(VECTOR *v, double mass, double temp)
 * 	Возвращает вектор v с компонентами распределенными в полупространство (x > 0) с температурой temp и массой mass.
 */
extern void get_vx_half_maxwell(VECTOR *v, double mass, double temp);

/*!
 * \fn void move_particle_to_buf(PARTICLE *p, PARTICLE *b, VECTOR *v, double new_weight)
 * 	Увеличивает вес частицы b на new_weight и с вероятностью new_weight / b->weight заносит в b частицу с массой, диаметром и скоростью p->mass, p->diam и v соответственно.
 */
extern void move_particle_to_buf(PARTICLE *p, PARTICLE *b, VECTOR *v, double new_weight);
/*!
 * \fn double gam(double x)
 * Возвращает гамма функцию от \param x относительно 0
 */
extern double gam(double x);

/*!
 * \fn double lbs(double XIA, double XIB)
 * Возвращает распределение внутренней энергии 
 */
extern double lbs(double XIA, double XIB);

#endif
