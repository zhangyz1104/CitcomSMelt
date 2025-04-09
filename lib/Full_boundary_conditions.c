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
#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

#include "lith_age.h"
#ifdef USE_GGRD
#include "ggrd_handling.h"
#endif
/* ========================================== */

void horizontal_bc(struct All_variables *, float *[], int, int, float,
                   unsigned int, char, int, int);
void assign_internal_bc(struct All_variables *);
void velocity_apply_periodic_bcs();
void temperature_apply_periodic_bcs();
void read_temperature_boundary_from_file(struct All_variables *);
void read_velocity_boundary_from_file(struct All_variables *);

/* ========================================== */

void full_velocity_boundary_conditions(E) struct All_variables *E;
{
  void velocity_imp_vert_bc();
  void velocity_apply_periodicapply_periodic_bcs();

  void apply_side_sbc();

  int j, noz, lv, k, node;

  for (lv = E->mesh.gridmax; lv >= E->mesh.gridmin; lv--)
    for (j = 1; j <= E->sphere.caps_per_proc; j++) {
      noz = E->mesh.NOZ[lv];
      if (E->mesh.topvbc != 1) { /* free slip top */
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 1, 0.0, VBX, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 3, 0.0, VBZ, 1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 2, 0.0, VBY, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 1, E->control.VBXtopval, SBX,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 3, 0.0, SBZ, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 2, E->control.VBYtopval, SBY,
                      1, lv, j);
#ifdef USE_GGRD
        /* Ggrd traction control */
        if ((lv == E->mesh.gridmax) && E->control.ggrd.vtop_control)
          ggrd_read_vtop_from_file(E, TRUE);
#endif
      }
      if (E->mesh.botvbc != 1) { /* free slip bottom */
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 1, 0.0, VBX, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 3, 0.0, VBZ, 1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 2, 0.0, VBY, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 1, E->control.VBXbotval, SBX,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 3, 0.0, SBZ, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 2, E->control.VBYbotval, SBY,
                      1, lv, j);
      }

      if (E->mesh.topvbc == 1) { /* velocity/no slip BC */
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 1, E->control.VBXtopval, VBX,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 3, 0.0, VBZ, 1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 2, E->control.VBYtopval, VBY,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 1, 0.0, SBX, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 3, 0.0, SBZ, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, noz, 2, 0.0, SBY, 0, lv, j);

#ifdef USE_GGRD
        /* Ggrd velocity (and age) control (only for the top processors */
        if ((lv == E->mesh.gridmax) && E->control.ggrd.vtop_control)
          ggrd_read_vtop_from_file(E, TRUE);
#endif

        if (E->control.vbcs_file) { /* this should either only be called
                                       once, or the input routines need
                                       to be told what to do for each
                                       multigrid level and cap. it might
                                       be easiest to call only once and
                                       have routines deal with multigrid
                                    */
          if ((lv == E->mesh.gridmin) && (j == E->sphere.caps_per_proc))
            read_velocity_boundary_from_file(E);
        }
      }

      if (E->mesh.botvbc == 1) { /* velocity bottom BC */
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 1, E->control.VBXbotval, VBX,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 3, 0.0, VBZ, 1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 2, E->control.VBYbotval, VBY,
                      1, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 1, 0.0, SBX, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 3, 0.0, SBZ, 0, lv, j);
        horizontal_bc(E, E->sphere.cap[j].VB, 1, 2, 0.0, SBY, 0, lv, j);
      }
    } /* end for j and lv */

  if (E->control.side_sbcs)
    apply_side_sbc(E);

  /* if(E->control.verbose) { */
  /*  for (j=1;j<=E->sphere.caps_per_proc;j++) */
  /*    for (node=1;node<=E->lmesh.nno;node++) */
  /*       fprintf(E->fp_out,"m=%d VB== %d %g %g %g flag %u %u
   * %u\n",j,node,E->sphere.cap[j].VB[1][node],E->sphere.cap[j].VB[2][node],E->sphere.cap[j].VB[3][node],E->node[j][node]&VBX,E->node[j][node]&VBY,E->node[j][node]&VBZ);
   */
  /*  fflush(E->fp_out); */
  /* } */

  /* If any imposed internal velocity structure it goes here */

  /*
    apply stress or velocity boundary conditions, read from file
    settings are to be implemented in those routines (will only do
    anything at present, if E->mesh.toplayerbc != 0
  */
  assign_internal_bc(E);

  return;
}

/* ========================================== */

void full_temperature_boundary_conditions(E) struct All_variables *E;
{
  void temperatures_conform_bcs();
  void temperature_imposed_vert_bcs();
  int j, lev, noz;
  if (E->mesh.toptbc_pole)
    if (!E->mesh.toptbc)
      myerror(E, "top pole temperature BC requires regular temperature BC");

  lev = E->mesh.levmax;
  for (j = 1; j <= E->sphere.caps_per_proc; j++) {
    noz = E->mesh.noz;
    if (E->mesh.toptbc == 1) {
      if (E->mesh.toptbc_pole) /* latitude dependency */
        horizontal_bc_lat_dep(E, E->sphere.cap[j].TB, noz, 3, TBZ, 1, lev, j);
      else
        horizontal_bc(E, E->sphere.cap[j].TB, noz, 3, E->control.TBCtopval, TBZ,
                      1, lev, j);
      horizontal_bc(E, E->sphere.cap[j].TB, noz, 3, E->control.TBCtopval, FBZ,
                    0, lev, j);
      if (E->control.tbcs_file)
        read_temperature_boundary_from_file(E);
    } else {
      horizontal_bc(E, E->sphere.cap[j].TB, noz, 3, E->control.TBCtopval, TBZ,
                    0, lev, j);
      horizontal_bc(E, E->sphere.cap[j].TB, noz, 3, E->control.TBCtopval, FBZ,
                    1, lev, j);
    }

    if (E->mesh.bottbc == 1) {
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, TBZ, 1,
                    lev, j);
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, FBZ, 0,
                    lev, j);
    } else {
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, TBZ, 0,
                    lev, j);
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, FBZ, 1,
                    lev, j);
    }
    if (E->control.lith_age_time == 1) {
      /* set the regions in which to use lithosphere files to determine
         temperature note that this is called if the lithosphere age in inputted
         every time step OR it is only maintained in the boundary regions */
      lith_age_temperature_bound_adj(E, lev);
    }
  } /* end for j */

  temperatures_conform_bcs(E);
  E->temperatures_conform_bcs = temperatures_conform_bcs;

  return;
}
/*  =========================================================  */

void horizontal_bc(struct All_variables *E, float *BC[], int ROW, int dirn,
                   float value, unsigned int mask, char onoff, int level,
                   int m) {
  int i, j, node, rowl;

  /* safety feature */
  if (dirn > E->mesh.nsd)
    return;

  if (ROW == 1)
    rowl = 1;
  else
    rowl = E->lmesh.NOZ[level];

  if (((ROW == 1) && (E->parallel.me_loc[3] == 0)) ||
      ((ROW == E->mesh.NOZ[level]) &&
       (E->parallel.me_loc[3] == E->parallel.nprocz - 1))) {

    /* turn bc marker to zero */
    if (onoff == 0) {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++) {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] +
                 (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][m][node] = E->NODE[level][m][node] & (~mask);
        } /* end for loop i & j */
    }

    /* turn bc marker to one */
    else {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++) {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] +
                 (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][m][node] = E->NODE[level][m][node] | (mask);
          if (level == E->mesh.levmax) /* NB */
            BC[dirn][node] = value;
        } /* end for loop i & j */
    }

  } /* end for if ROW */

  return;
}

/* latitudinally dependent horizontal BC */

void horizontal_bc_lat_dep(struct All_variables *E, float *BC[], int ROW,
                           int dirn, unsigned int mask, char onoff, int level,
                           int m) {
  int i, j, node, rowl;

  /* safety feature */
  if (dirn > E->mesh.nsd)
    return;

  if (ROW == 1)
    rowl = 1;
  else
    rowl = E->lmesh.NOZ[level];

  if (((ROW == 1) && (E->parallel.me_loc[3] == 0)) ||
      ((ROW == E->mesh.NOZ[level]) &&
       (E->parallel.me_loc[3] == E->parallel.nprocz - 1))) {

    /* turn bc marker to zero */
    if (onoff == 0) {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++) {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] +
                 (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][m][node] = E->NODE[level][m][node] & (~mask);
        } /* end for loop i & j */
    }

    /* turn bc marker to one */
    else {
      for (j = 1; j <= E->lmesh.NOY[level]; j++)
        for (i = 1; i <= E->lmesh.NOX[level]; i++) {
          node = rowl + (i - 1) * E->lmesh.NOZ[level] +
                 (j - 1) * E->lmesh.NOX[level] * E->lmesh.NOZ[level];
          E->NODE[level][m][node] = E->NODE[level][m][node] | (mask);
          if (level == E->mesh.levmax) /* NB */
            BC[dirn][node] = lat_dep_temp(E, node);
        } /* end for loop i & j */
    }

  } /* end for if ROW */

  return;
}

/* temp as f(latitude) */
float lat_dep_temp(struct All_variables *E, int node) {
  float xp[3], fac;
  xyz2rtp(E->x[1][1][node], E->x[1][2][node], E->x[1][3][node], xp);
  fac = sin(xp[1]); /* sin(theta), 0 at pole, 1  equator */
  return E->control.TBCtop_pole +
         fac * (E->control.TBCtopval - E->control.TBCtop_pole);
}

void velocity_apply_periodic_bcs(E) struct All_variables *E;
{
  fprintf(E->fp, "Periodic boundary conditions\n");

  return;
}

void temperature_apply_periodic_bcs(E) struct All_variables *E;
{
  fprintf(E->fp, "Periodic temperature boundary conditions\n");

  return;
}

/* update cmb temperature by average heat flux on cmb */
void update_cmb_T(struct All_variables *E) {
  double rhom, rhoc, Cm, Cc, rc, deltat, gradT, deltaT, gradT_melt;
  int j, lev;

  lev = E->mesh.levmax;
  rhom = E->data.density;
  rhoc = E->data.density_below;
  Cm = E->data.Cp;
  Cc = E->mush.dynamics.Cpc;
  rc = E->sphere.ri;
  gradT = E->sphere.cap[1].heat_flux;
  gradT_melt = E->sphere.cap[1].melt_flux;
  deltat = E->advection.timestep;

  // set segregation heatflux to zero
  // gradT_melt = 0;

  if (E->parallel.me == 0)
    fprintf(stderr, "TBCbotval = %.4g, gradT = %.3g, gradT_melt = %.3g\n",
            E->control.TBCbotval, gradT, gradT_melt);

  /*see nanzhang etal 2013 for reference, undimentionalize all components in eq 12 */
  deltaT = 0;
  if (E->mush.dynamics.convection_developed > 0)
    deltaT = -3 * rhom * Cm / (rhoc * Cc) / rc * (gradT + gradT_melt) * deltat;
  E->control.TBCbotval += (float) deltaT;

  // ZYZ debug
  // E->control.TBCbotval = E->Have.T[2];

  /*change and conform boundary condition with updated cmb temperature*/
  for (j = 1; j <= E->sphere.caps_per_proc; j++)
    if (E->mesh.bottbc == 1) {
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, TBZ, 1, lev, j);
      horizontal_bc(E, E->sphere.cap[j].TB, 1, 3, E->control.TBCbotval, FBZ, 0, lev, j);
    }
  temperatures_conform_bcs(E);
}