#ifndef INITIAL_STATES_H
#define INITIAL_STATES_H

int init();

/*! Инициализация камеры 2 потоками. Частицы с номером 0 в ячейке движутся вниз, 1 - вверх.
*/
int init_2flows();

void read_data();

int create_camera();

int init_from_file(char *file_name);

#endif
