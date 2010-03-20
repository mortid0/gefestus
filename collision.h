#ifndef COLLISION_H_
#define COLLISION_H_

extern void elastic_collision(int r);

/*!
 * \param BN is coefficient of KT in the mean of translational relative energy
 */
#define BN 2.0

extern void inelastic_collision(int krp);
#endif
