#ifndef PICTURES_H
#define PICTURES_H

/* О©╫О©╫О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫. О©╫О©╫О©╫О©╫ О©╫О©╫ О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫ - О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫ О©╫О©╫О©╫О©╫О©╫О©╫О©╫О©╫ 0 */

int SHOW_TEMP_FREQ;
int SHOW_CONC_FREQ;
int SHOW_MASS_FREQ;
#define SHOW_NALIP_FREQ 0
int SHOW_IMPULSE_FREQ;

/*! \fn void show_field(int num, double (*cell_func)(int, int, int), char *func_dir, int startx, int endx, int starty, int endy);
 * Функция печатает поле распределения заданной функции газа(тепературы, концентрации и .д.)
 * \param num - номер расчетного шага
 * \param cell_func(nc,x,y) - функция газа в заданной точке камере(x,y) для определенного компонента(nc)
 * \param func_dir - название каталога, в который будет сохраняться снимок.
 * \param startx - номер ячейки по горизонтали, с которой начинается съем параметра
 * \param endx - номер ячейки по горизонтали, на которой заканчивается съем параметра
 * \param starty - номер ячейки по вертикали, с которой начинается съем параметра
 * \param endy - номер ячейки по вертикали, на которой заканчивается съем параметра
 */
extern void show_field(int num, double (*cell_func)(int, int, int), char *func_dir, int startx, int endx, int starty, int endy);

/*! Выводит поле температур для 2 потоков. 1 столбец - частицы вниз, 2 - вверх.
 */
extern void show_temp_2flows(int num);

/*!Выводит поле концентраций для 2 потоков. 1 столбец - частицы вниз, 2 - вверх.
 */
extern void show_conc_2flows(int num);

/*!Выводит поле vy для 2 потоков.
 */
extern void show_vy_2flows(int num);

extern void show_vector(int num);

//extern void show_mass(int num);

//extern void show_nalip(int num);

/*!
 * Снимает необходимые параметры и сохраняет их на диске.
 * \param step - номер текущего расчетного шага.
*/
extern void print_pictures(int step);

#endif
