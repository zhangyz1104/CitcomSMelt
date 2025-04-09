/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/*  Here are the routines which process the results of each buoyancy solution,
   and call any relevant output routines. Much of the information has probably
   been output along with the velocity field. (So the velocity vectors and other
   data are fully in sync). However, heat fluxes and temperature averages are
   calculated here (even when they get output the next time around the velocity
   solver);
    */

#include "element_definitions.h"
#include "global_defs.h"
#include <math.h> /* for sqrt */

void parallel_process_termination(void);

void post_processing(struct All_variables *E) { return; }

/* ===================
    Surface heat flux
   =================== */

void heat_flux(E)struct All_variables *E;
{
  int m, e, el, i, j, node, lnode, nz;
  float *flux[NCS], *SU[NCS], *RU[NCS];
  float VV[4][9], u[9], T[9], dTdz[9], rho[9], uT;
  float *sum_h;

  void velo_from_element();
  void sum_across_surface();
  void return_horiz_ave();
  void return_horiz_ave_f();

  const int dims = E->mesh.nsd, dofs = E->mesh.dof;
  const int vpts = vpoints[dims];
  const int ppts = ppoints[dims];
  const int ends = enodes[dims];
  const int nno = E->lmesh.nno;
  const int lev = E->mesh.levmax;
  const int sphere_key = 1;

  sum_h = (float *) malloc((5) * sizeof(float));
  for (i = 0; i <= 4; i++)
    sum_h[i] = 0.0;

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {

    flux[m] = (float *) malloc((1 + nno) * sizeof(float));

    for (i = 1; i <= nno; i++) {
      flux[m][i] = 0.0;
    }

    for (e = 1; e <= E->lmesh.nel; e++) {

      velo_from_element(E, VV, m, e, sphere_key);

      for (i = 1; i <= vpts; i++) {
        u[i] = 0.0;
        T[i] = 0.0;
        dTdz[i] = 0.0;
        rho[i] = 0.0;
        for (j = 1; j <= ends; j++) {
          nz = ((E->ien[m][e].node[j] - 1) % E->lmesh.noz) + 1;
          rho[i] += E->refstate.rho[nz] * E->N.vpt[GNVINDEX(j, i)];
          u[i] += VV[3][j] * E->N.vpt[GNVINDEX(j, i)];
          T[i] += E->T[m][E->ien[m][e].node[j]] * E->N.vpt[GNVINDEX(j, i)];
          dTdz[i] += -E->T[m][E->ien[m][e].node[j]] * E->gNX[m][e].vpt[GNVXINDEX(2, j, i)];
        }
      }

      uT = 0.0;
      /* XXX: missing unit conversion, heat capacity and thermal conductivity */
      for (i = 1; i <= vpts; i++)
        uT += rho[i] * u[i] * T[i] * E->gDA[m][e].vpt[i] + dTdz[i] * E->gDA[m][e].vpt[i];

      uT /= E->eco[m][e].area;
      for (j = 1; j <= ends; j++)
        flux[m][E->ien[m][e].node[j]] += uT * E->TWW[lev][m][e].node[j];

    } /* end of e */
  }   /* end of m */

  (E->exchange_node_f)(E, flux, lev);

  for (m = 1; m <= E->sphere.caps_per_proc; m++)
    for (i = 1; i <= nno; i++)
      flux[m][i] *= E->MASS[lev][m][i];

  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      for (i = 1; i <= E->lmesh.nsf; i++)
        E->slice.shflux[m][i] = flux[m][E->surf_node[m][i]];

  if (E->parallel.me_loc[3] == 0)
    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      for (i = 1; i <= E->lmesh.nsf; i++)
        E->slice.bhflux[m][i] = flux[m][E->surf_node[m][i] - E->lmesh.noz + 1];

  for (m = 1; m <= E->sphere.caps_per_proc; m++)
    for (e = 1; e <= E->lmesh.snel; e++) {
      uT = (E->slice.shflux[m][E->sien[m][e].node[1]] +
          E->slice.shflux[m][E->sien[m][e].node[2]] +
          E->slice.shflux[m][E->sien[m][e].node[3]] +
          E->slice.shflux[m][E->sien[m][e].node[4]]) *
          0.25;
      el = e * E->lmesh.elz;
      sum_h[0] += uT * E->eco[m][el].area;
      sum_h[1] += E->eco[m][el].area;

      uT = (E->slice.bhflux[m][E->sien[m][e].node[1]] +
          E->slice.bhflux[m][E->sien[m][e].node[2]] +
          E->slice.bhflux[m][E->sien[m][e].node[3]] +
          E->slice.bhflux[m][E->sien[m][e].node[4]]) *
          0.25;
      el = (e - 1) * E->lmesh.elz + 1;
      sum_h[2] += uT * E->eco[m][el].area;
      sum_h[3] += E->eco[m][el].area;
    }

  sum_across_surface(E, sum_h, 4);

  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
    sum_h[0] = sum_h[0] / sum_h[1];

    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      E->sphere.cap[m].heat_flux = sum_h[0];  // ZYZ-220218

    if (E->parallel.me == E->parallel.nprocz - 1)
      fprintf(stderr, "surface heat flux= %f\n", sum_h[0]);
  }

  if (E->parallel.me_loc[3] == 0) {
    sum_h[2] = sum_h[2] / sum_h[3];

    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      E->sphere.cap[m].heat_flux = sum_h[2];  // ZYZ-220218

    if (E->parallel.me == 0)
      fprintf(stderr, "bottom heat flux= %f\n", sum_h[2]);
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++)
    free((void *) flux[m]);

  free((void *) sum_h);

  return;
}

/*
  compute horizontal average of temperature, composition and rms velocity
  add viscosity and melt fraction (ZYZ-220221)
*/
void compute_horiz_avg(struct All_variables *E) {
  void return_horiz_ave_f();

  int m, n, i;
  float *S1[NCS], *S2[NCS], *S3[NCS], *S4[NCS];
  float *S5[NCS], *S6[NCS], *S7[NCS], *S8[NCS], *S9[NCS], *S10[NCS], *S11[NCS];

  const double L = E->mush.consts.St0;

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    S1[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S2[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S3[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S4[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S5[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S6[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S7[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S8[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S9[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S10[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
    S11[m] = (float *) malloc((E->lmesh.nno + 1) * sizeof(float));
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    for (i = 1; i <= E->lmesh.nno; i++) {
      S1[m][i] = E->T[m][i];
      S2[m][i] = E->sphere.cap[m].V[1][i] * E->sphere.cap[m].V[1][i] +
          E->sphere.cap[m].V[2][i] * E->sphere.cap[m].V[2][i];
      S3[m][i] = E->sphere.cap[m].V[3][i];
      S4[m][i] = E->VI[E->mesh.levmax][m][i];
      S5[m][i] = E->mush.comp.phi[i];
      S6[m][i] = E->mush.dynamics.q[i];
      S7[m][i] = E->mush.comp.c0[E->mush.consts.ncomp][i];
      S8[m][i] = E->mush.comp.cl[E->mush.consts.ncomp][i];
      S9[m][i] = E->mush.dynamics.E[i];
//      S10[m][i] = E->mush.dynamics.M[i];
//      S11[m][i] = E->mush.comp.drho2[i];
      S10[m][i] = S3[m][i] * (S1[m][i] + L * S5[m][i]);
      S11[m][i] = S6[m][i] * (S1[m][i] + L);
    }
  }

  return_horiz_ave_f(E, S1, E->Have.T);
  return_horiz_ave_f(E, S2, E->Have.V[1]);
  return_horiz_ave_f(E, S3, E->Have.V[2]);
  return_horiz_ave_f(E, S4, E->Have.vis);
  return_horiz_ave_f(E, S5, E->Have.phi);
  return_horiz_ave_f(E, S6, E->Have.q);
  return_horiz_ave_f(E, S7, E->Have.c0);
  return_horiz_ave_f(E, S8, E->Have.cl);
  return_horiz_ave_f(E, S9, E->Have.E);
  return_horiz_ave_f(E, S10, E->Have.M);
  return_horiz_ave_f(E, S11, E->Have.drho2);

  if (E->composition.on) {
    for (n = 0; n < E->composition.ncomp; n++) {
      for (m = 1; m <= E->sphere.caps_per_proc; m++) {
        for (i = 1; i <= E->lmesh.nno; i++)
          S1[m][i] = E->composition.comp_node[m][n][i];
      }
      return_horiz_ave_f(E, S1, E->Have.C[n]);
    }
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    free((void *) S1[m]);
    free((void *) S2[m]);
    free((void *) S3[m]);
    free((void *) S4[m]);
    free((void *) S5[m]);
    free((void *) S6[m]);
    free((void *) S7[m]);
    free((void *) S8[m]);
    free((void *) S9[m]);
    free((void *) S10[m]);
    free((void *) S11[m]);
  }

  for (i = 1; i <= E->lmesh.noz; i++)
    E->Have.V[1][i] = sqrt(E->Have.V[1][i]);

//  if (E->mush.mush_version)
//    composition_conservation(E);  // mush-related
}
