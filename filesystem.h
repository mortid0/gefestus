#ifndef FILESYSTEM_H
#define FILESYSTEM_H


/*!
 * Создает структуру каталогов для работы комплекса.
 */
extern void make_dirs();

/*!
 * Сохраняет текущее состояние камеры и переменных на расчетном шаге \param step в файл на диск.
 */
extern void save_camera(int step);

/*!
 * Загружает состояние камеры и переменных с расчетного шага \param step
 */
extern void load_camera(int step);
#endif
