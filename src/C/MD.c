//  Performance Programming
//  Coursework
//  B063622


/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <math.h>
#include "coord.h"

void visc_force(int N, double* f, double* visc, double* vel);
void add_norm(int N, double* r, double* delta);
double force(double W, double delta, double r);
void wind_force(int N, double* f, double* visc, double vel);


void evolve(int count, double dt)
{
	int  step;
	int i, j, k, l;
	int collided;

	/*
	 * Loop over timesteps.
	 */
	for (step = 1; step <= count; step++)
	{
		printf("timestep %d\n", step);
		printf("collisions %d\n", collisions);

		/* set the viscosity term in the force calculation */
		/* add the wind term in the force calculation */
		double *pf    = f[0];
		double *pvel  = vel[0];
		double *pwind = wind;
		
		for (j = 0; j < Ndim; j++)
		{
			double *pvi = visc;
			
			for (i = 0; i < Nbody; i++)
			{
				*pf = -(*pvi) * (*pvel + *pwind);
				
				++pvi;
				++pf;
				++pvel;
			}
			
			++pwind;
		}

		/* calculate distance from central mass */
		for (k = 0; k < Nbody; k++)
		{
			r[k] = 0.0;

			for (i = 0; i < Ndim; i++)
			{
				r[k] += (pos[i][k] * pos[i][k]);
			}
			
			r[k] = sqrt(r[k]);
		}

		/* calculate central force */
		for (l = 0; l < Ndim; l++)
		{
			for (i = 0; i < Nbody; i++)
			{
				f[l][i] = f[l][i] - force(G * mass[i] * M_central, pos[l][i], r[i]);
			}
		}

		/* calculate pairwise separation of particles */
		for (l = 0; l < Ndim; l++)
		{
			k = 0;
			for (i = 0; i < Nbody; i++)
			{
				for (j = i + 1; j < Nbody; j++)
				{
					delta_pos[l][k] = pos[l][i] - pos[l][j];
					k = k + 1;
				}
			}
		}

		/* calculate norm of seperation vector */
		for (k = 0; k < Npair; k++)
		{
			delta_r[k] = 0.0;
		}
		for (i = 0; i < Ndim; i++)
		{
			add_norm(Npair, delta_r, delta_pos[i]);
		}
		for (k = 0; k < Npair; k++)
		{
			delta_r[k] = sqrt(delta_r[k]);
		}

		/*
		 * add pairwise forces.
		 */
		int local_collisions = 0;
		for (l = 0; l < Ndim; l++)
		{
			k = 0;
			double *pfli = f[l];			
			for (i = 0; i < Nbody; i++)
			{
				double *pflj = &f[l][i+1];
				for (j = i + 1; j < Nbody; j++)
				{
					/*  
					 * flip force if close in
					 */
					double multiplier = 1.0;
					if (! (delta_r[k] >= Size))
					{
						multiplier = -1.0;
						local_collisions++;
					}
						
					double gforce = multiplier * G * mass[i] * mass[j] * delta_pos[l][k] / pow(delta_r[k], 3.0);
					
					*pfli -= gforce;					
					*pflj += gforce;
					
					++k;
					++pflj;
				}
				
				++pfli;
			}
		}

		//  this is a kludge: the modified code results in 3x as many collisions so we correct this here
		collisions += local_collisions / 3;
		
		/* update positions */
		double *ppos = pos[0];
		pvel         = vel[0];
		
		for (j = 0; j < Ndim; j++)
		{
			for (i = 0; i < Nbody; i++)
			{
				*ppos += dt * *pvel;
				++ppos;
				++pvel;
			}
		}

		/* update velocities */
		pvel = vel[0];
		pf  = f[0];
		
		for (j = 0; j < Ndim; j++)
		{
			double *pmass = mass;
			for (i = 0; i < Nbody; i++)
			{
				*pvel += dt * (*pf / *pmass);
				++pvel;
				++pf;
				++pmass;
			}
		}
	}
}




