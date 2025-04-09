#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "global_defs.h"
#include "element_definitions.h"

static void Tm(struct All_variables *E);
static double K(struct All_variables *E, double T, int i, int n);
static double residual(struct All_variables *E, double T, double phi, int i);
static double Newton_iteration_T(struct All_variables *E, double phi, int i);
static double Newton_iteration_phi(struct All_variables *E, double T, int i);
static void solid_flow(struct All_variables *E, int el, double dc0[9], int n);

/* compute melt fraction, melt composition, density and Stefan Number */
void thermodynamic_calculation(struct All_variables *E) {
  int i, n;
  double T, phi, Kn, dT;
  double cl, cs, rhol, rhos;

  const int m = 1;
  const int nno = E->lmesh.nno;
  const int nc = E->mush.consts.ncomp;
  const double St0 = E->mush.consts.St0 * E->data.ref_temperature;
  const double rho0 = E->data.density;
  const double drho0 = E->mush.consts.density_contrast;
  const double phic = 0.5;

  dT = 0.01;  // St_max=100

  // init melting point
  if (E->monitor.solution_cycles == 0)
    Tm(E);

  // compute solidus Tsol(P, C) and liquidus Tliq(P, C)
  for (i = 1; i <= nno; i++) {
    E->mush.comp.Tsol[i] = Newton_iteration_T(E, 0, i);
    E->mush.comp.Tliq[i] = Newton_iteration_T(E, 1, i);
  }

  // compute melt fraction phi(T, P, C) and composition c^{s, l}(T, P, C)
  for (i = 1; i <= nno; i++) {
    // melt fraction
    T = E->T[m][i] * E->data.ref_temperature;
    phi = max(0, min(1, Newton_iteration_phi(E, T, i)));
//    if (phi > phic) {
//      T = Newton_iteration_T(E, phic, i);
//      E->T[m][i] = T / E->data.ref_temperature;
//      phi = phic;
//    }
    E->mush.comp.phi[i] = phi;

    // mass concentration
    for (n = 1; n <= nc; n++) {
      Kn = K(E, T, i, n);
      cl = E->mush.comp.c0[n][i] / (phi + (1 - phi) * Kn);
      E->mush.comp.cl[n][i] = max(0, min(1, cl));
      cs = E->mush.comp.c0[n][i] / (phi / Kn + (1 - phi));
      E->mush.comp.cs[n][i] = max(0, min(1, cs));
    }

    // density
    rhol = rhos = 0;
    for (n = 1; n <= nc; n++) {
      rhol += E->mush.comp.cl[n][i] / E->mush.consts.melt_density[n];
      rhos += E->mush.comp.cs[n][i] / E->mush.consts.solid_density[n];
    }
    rhol = 1 / rhol;
    rhos = 1 / rhos;
    if (phi == 0)
      rhol = E->mush.consts.melt_density[2];
    else if (phi == 1)
      rhos = E->mush.consts.solid_density[1];
//    if (E->mush.comp.c0[1][i] == 1) {
//      rhol = E->mush.consts.melt_density[1];
//      rhos = E->mush.consts.solid_density[1];
//    } else if (E->mush.comp.c0[2][i] == 1) {
//      rhol = E->mush.consts.melt_density[2];
//      rhos = E->mush.consts.solid_density[2];
//    }

    E->mush.comp.drho[i] = (rhos - rhol) / drho0;
    E->mush.comp.drho2[i] = -(phi * rhol + (1 - phi) * rhos - rho0) / drho0;
    if (E->mush.comp.constant_density) {
      E->mush.comp.drho[i] = 1;
      E->mush.comp.drho2[i] = phi;
    }

    // Stefan number
    T -= dT;
    if (T <= E->mush.comp.Tsol[i])
      phi = 0;
    else
      phi = max(0, min(1, Newton_iteration_phi(E, T, i)));
    E->mush.comp.St[i] = St0 * (E->mush.comp.phi[i] - phi) / dT;

    if (isnan(E->mush.comp.drho2[i]) || isinf(E->mush.comp.drho2[i]) || isnan(E->T[m][i]) || isinf(E->T[m][i])) {
      fprintf(stderr, "\n\n\n--------debug information2---------\n");
      fprintf(stderr,
              "T=%.4e, phi=%.4e, Tsol=%.4e, Tliq=%.4e, c^ol=%.4e, c^an=%.4e, drho=%.4e, rhol=%.4e, rhos=%.4e, St=%.4e\n",
              T,
              E->mush.comp.phi[i],
              E->mush.comp.Tsol[i],
              E->mush.comp.Tliq[i],
              E->mush.comp.c0[1][i],
              E->mush.comp.c0[2][i],
              E->mush.comp.drho2[i],
              rhol,
              rhos,
              E->mush.comp.St[i]);
    }
  }
}

void update_bulk_composition(struct All_variables *E) {
  int i, j, n, el, node, isurf;
  double c0, dqcldr, dqdr, dV, dNadr, PG;
  double *dc0_darcy[NCS], *dc0_solid[NCS], dc0_dike, dc0[9];

  const int m = 1;
  const int nel = E->lmesh.nel;
  const int nno = E->lmesh.nno;
  const int noz = E->lmesh.noz;
  const int lev = E->mesh.levmax;
  const int nc = E->mush.consts.ncomp;
  const float dt = E->advection.timestep;
  const double dr_surf = E->eco[1][noz - 1].size[3];  // surface grid radius

  dc0_darcy[m] = (double *) malloc((nno + 1) * sizeof(double));
  dc0_solid[m] = (double *) malloc((nno + 1) * sizeof(double));

  // update surface composition
  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    for (isurf = 1; isurf <= E->lmesh.nsf; isurf++) {
      i = isurf * noz;
      if (E->mush.dynamics.M[i])
        for (n = 1; n <= nc; n++)
          E->mush.comp.c0[n][i] += E->mush.dynamics.Mn[n][i] * dt / dr_surf;
    }

  // update internal composition:
  // dc0/dt = -vs * grad(c0) -(dqcl/dr - c0 * dq/dr) -E(cl - c0)
  for (n = 1; n <= nc; n++) {
    for (i = 1; i <= nno; i++)
      dc0_darcy[m][i] = dc0_solid[m][i] = 0;

    // melt segregation: + solid flow
    for (el = 1; el <= nel; el++) {
      // melt segregation: dqcl/dr - c0 * dq/dr
      for (i = 1; i <= 8; i++) {
        c0 = dqdr = dqcldr = 0;
        for (j = 1; j <= 8; j++) {
          node = E->ien[m][el].node[j];
          dNadr = E->gNX[m][el].vpt[GNVXINDEX(2, j, i)];
          c0 += E->mush.comp.c0[n][node] * E->N.vpt[GNVINDEX(j, i)];
          dqcldr += E->mush.dynamics.q[node] * E->mush.comp.cl[n][node] * dNadr;
          dqdr += E->mush.dynamics.q[node] * dNadr;
        }
        dc0[i] = dqcldr - c0 * dqdr;
      }

      for (j = 1; j <= 8; j++)
        for (i = 1; i <= 8; i++) {
          dV = E->gDA[m][el].vpt[i];
          PG = E->N.vpt[GNVINDEX(j, i)] + E->eco[m][el].size[3] * E->gNX[m][el].vpt[GNVXINDEX(2, j, i)] * 0.5;
          dc0_darcy[m][E->ien[m][el].node[j]] += dc0[i] * PG * dV;
        }

      // solid flow: vs * grad(c0) (SUPG scheme)
      solid_flow(E, el, dc0, n);
      for (j = 1; j <= 8; j++)
        dc0_solid[m][E->ien[m][el].node[j]] += dc0[j];
    }

    (E->exchange_node_d)(E, dc0_darcy, lev);
    (E->exchange_node_d)(E, dc0_solid, lev);

    for (i = 1; i <= nno; i++) {
      dc0_solid[m][i] *= E->Mass[m][i];
      dc0_darcy[m][i] *= E->Mass[m][i];
      if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
        dc0_dike = E->mush.dynamics.E[i] * (E->mush.comp.cl[n][i] - E->mush.comp.c0[n][i]);
      else
        dc0_dike = 0;

      // no-flux lower boundary
      if ((E->mush.dynamics.solve_CMB_comp) && (E->node[m][i] & TBZ) && (E->parallel.me_loc[3] == 0))
        dc0_darcy[m][i] = E->mush.dynamics.q[i] * (E->mush.comp.cl[n][i] - E->mush.comp.c0[n][i]) / E->eco[1][1].size[3];

      // update dc0
      E->mush.comp.c0[n][i] -= (dc0_solid[m][i] + dc0_darcy[m][i] + dc0_dike) * dt;
    }
  }

  // norm
  for (i = 1; i <= nno; i++) {
    c0 = 0;
    for (n = 1; n <= nc; n++) {
      E->mush.comp.c0[n][i] = max(E->mush.comp.c0[n][i], 0);
      c0 += E->mush.comp.c0[n][i];
    }
    for (n = 1; n <= nc; n++)
      E->mush.comp.c0[n][i] /= c0;
  }

  free(dc0_darcy[m]);
  free(dc0_solid[m]);
}

void update_bulk_composition_old(struct All_variables *E) {
  int i, j, n, el, node, isurf, diff_nodes;
  double c0, dqcldr, dqdr, dV, dNadr, PG;
  double dr_dike, c_dike;
  double *dc0_darcy[NCS], *dc0_solid[NCS], dc0_dike, dc0[9];

  const int m = 1;
  const int nel = E->lmesh.nel;
  const int nno = E->lmesh.nno;
  const int noz = E->lmesh.noz;
  const int lev = E->mesh.levmax;
  const int nc = E->mush.consts.ncomp;
  const float dt = E->advection.timestep;
  const double dr_surf = E->eco[1][32].size[3];  // surface grid radius

  dc0_darcy[m] = (double *) malloc((nno + 1) * sizeof(double));
  dc0_solid[m] = (double *) malloc((nno + 1) * sizeof(double));

  // update surface composition
  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    for (isurf = 1; isurf <= E->lmesh.nsf; isurf++) {
      i = isurf * noz;
      if (E->mush.dynamics.M[i]) {
//        dr_dike = E->mush.dynamics.M[i] * dt;
        for (n = 1; n <= nc; n++) {
//          c_dike = E->mush.dynamics.Mn[n][i] / E->mush.dynamics.M[i];
//          c0 = E->mush.comp.c0[n][i];
//          c0 += (c_dike - c0) * dr_dike / dr_surf;
//          E->mush.comp.c0[n][i] = E->mush.comp.c0[n][i - 1] = c0;
//          if (dr_dike > dr_surf)
//            fprintf(stderr, "(!!!!!) melt extraction error: %.4g\n", dr_dike / dr_surf);
          E->mush.comp.c0[n][i] += E->mush.dynamics.Mn[n][i] * dt / dr_surf;
        }
      }
    }

  // update internal composition:
  // dc0/dt = -vs * grad(c0) -(dqcl/dr - c0 * dq/dr) -E(cl - c0)
  for (n = 1; n <= nc; n++) {
    for (i = 1; i <= nno; i++)
      dc0_darcy[m][i] = dc0_solid[m][i] = 0;

    // melt segregation: + solid flow
    for (el = 1; el <= nel; el++) {
      // melt segregation: dqcl/dr - c0 * dq/dr
      dc0[0] = 0;
      if ((el % E->lmesh.elz == 1) && (E->parallel.me_loc[3] == 0))
        diff_nodes = 5;
      else
        diff_nodes = 1;
      for (i = 1; i <= 8; i++) {
        c0 = dqdr = dqcldr = 0;
        for (j = diff_nodes; j <= 8; j++) {
          node = E->ien[m][el].node[j];
          dNadr = E->gNX[m][el].vpt[GNVXINDEX(2, j, i)];
          c0 += E->mush.comp.c0[n][node] * E->N.vpt[GNVINDEX(j, i)];
          dqcldr += E->mush.dynamics.q[node] * E->mush.comp.cl[n][node] * dNadr;
          dqdr += E->mush.dynamics.q[node] * dNadr;
        }
//        dV = E->gDA[m][el].vpt[i] / E->eco[m][el].area;
//        dc0[0] += (dqcldr - c0 * dqdr) * dV;
        dc0[i] += dqcldr - c0 * dqdr;
      }
//      for (j = 5; j <= 8; j++) // upwind scheme: only consider upward nodes
//        dc0_darcy[m][E->ien[m][el].node[j]] += dc0[0] * E->TWW[lev][m][el].node[j];

      for (j = 1; j <= 8; j++)
        for (i = 1; i <= 8; i++) {
          dV = E->gDA[m][el].vpt[i];
          PG = E->N.vpt[GNVINDEX(j, i)] + E->eco[m][el].size[3] * E->gNX[m][el].vpt[GNVXINDEX(2, j, i)] * 0.5;
          dc0_darcy[m][E->ien[m][el].node[j]] += dc0[i] * PG * dV;
        }

      // solid flow: vs * grad(c0) (SUPG scheme)
      solid_flow(E, el, dc0, n);
      for (j = 1; j <= 8; j++)
        dc0_solid[m][E->ien[m][el].node[j]] += dc0[j];
    }

    (E->exchange_node_d)(E, dc0_darcy, lev);
    (E->exchange_node_d)(E, dc0_solid, lev);

    for (i = 1; i <= nno; i++) {
      // no-flux lower boundary
      if ((E->mush.dynamics.solve_CMB_comp) && (E->node[m][i] & TBZ) && (E->parallel.me_loc[3] == 0)) {
        dc0_darcy[m][i] = dc0_darcy[m][i + 1] * 0.5;
        dc0_darcy[m][i + 1] *= 0.5;
      }

      dc0_solid[m][i] *= E->Mass[m][i];
      dc0_darcy[m][i] *= E->Mass[m][i];
      if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
        dc0_dike = E->mush.dynamics.E[i] * (E->mush.comp.cl[n][i] - E->mush.comp.c0[n][i]);
      else
        dc0_dike = 0;
      // dc0_dike = E->mush.dynamics.E[i] * E->mush.comp.cl[n][i] - E->mush.dynamics.Mn[n][i];

      // update dc0
      E->mush.comp.c0[n][i] -= (dc0_solid[m][i] + dc0_darcy[m][i] + dc0_dike) * dt;
      // E->mush.comp.c0[n][i] -= (dc0_solid[m][i] + dc0_dike) * dt;

      // output debug
      if (((E->parallel.me == 1) && (i <= 33) && (i >= 23)) && (n == 2))
        fprintf(stderr, "top[%d]: c0=%6.4g darcy=%8.2g solid=%8.2g dike=%8.2g "
                        "T=%6.4g phi=%8.4g E=%8.2g V=%8.2g\n",
                i, E->mush.comp.c0[n][i],
                -dc0_darcy[m][i] * dt, -dc0_solid[m][i] * dt, -dc0_dike * dt,
                E->T[m][i], E->mush.comp.phi[i],
                E->mush.dynamics.E[i] * dt, E->sphere.cap[m].V[3][i]);

      if ((E->parallel.me == 0) && (i <= 5) && (n == 2))
        fprintf(stderr, "    bottom[%d]: c0=%6.4g darcy=%8.2g solid=%8.2g "
                        "T=%4.2g phi=%8.2g\n",
                i, E->mush.comp.c0[n][i], -dc0_darcy[m][i] * dt, -dc0_solid[m][i] * dt,
                E->T[m][i], E->mush.comp.phi[i]);

      // output the velocity and segregation flux in one element for debugging the solid compaction
      if (E->mush.dynamics.solid_compaction && (E->parallel.me == 3) && (i == 48)) {
        fprintf(stderr, "[debug solid compaction]\n");
        fprintf(stderr, "x1=%.4g, x2=%.4g, x3=%.4g\n",
                E->eco[m][i].size[1], E->eco[m][i].size[2], E->eco[m][i].size[3]);
        for (j = 1; j <= 8; j++) {
          node = E->ien[m][i].node[j];
          fprintf(stderr, "%.4g, %.4g, %.4g, %.4g, ",
                  E->sphere.cap[m].V[1][node], E->sphere.cap[m].V[2][node],
                  E->sphere.cap[m].V[3][node], E->mush.dynamics.q[node]);
        }
        fprintf(stderr, "\n");
      }
    }
  }

  // norm
  for (i = 1; i <= nno; i++) {
    c0 = 0;
    for (n = 1; n <= nc; n++) {
      E->mush.comp.c0[n][i] = max(E->mush.comp.c0[n][i], 0);
      c0 += E->mush.comp.c0[n][i];
    }
    for (n = 1; n <= nc; n++)
      E->mush.comp.c0[n][i] /= c0;
  }

  free(dc0_darcy[m]);
  free(dc0_solid[m]);
}

void composition_conservation(struct All_variables *E) {
  int i;
  double dr, dV, V, V_all, dc0, c0, c0_all;
  double delta1, delta2;

  c0 = 0;
  V = 0;
  for (i = 2; i <= E->lmesh.noz; i++) {
    dr = E->sx[1][3][i] - E->sx[1][3][i - 1];
    dV = pow(E->sx[1][3][i] + E->sx[1][3][i - 1], 2) * dr;
    dc0 = (E->Have.c0[i] + E->Have.c0[i - 1]) * 0.5;
    V += dV;
    c0 += dc0 * dV; // only applied to two-component system
  }

  MPI_Allreduce(&c0, &c0_all, 1, MPI_DOUBLE, MPI_SUM, E->parallel.vertical_comm);
  MPI_Allreduce(&V, &V_all, 1, MPI_DOUBLE, MPI_SUM, E->parallel.vertical_comm);
  c0_all /= V_all;
  delta1 = 1 - c0_all - E->mush.consts.bulk_comp[1];
  delta2 = c0_all - E->mush.consts.bulk_comp[2];

  for (i = 1; i <= E->lmesh.nno; i++) {
    E->mush.comp.c0[1][i] = max(0, E->mush.comp.c0[1][i] - delta1);
    E->mush.comp.c0[2][i] = max(0, E->mush.comp.c0[2][i] - delta2);
  }

  if (E->parallel.me == 1)
    fprintf(stderr, "[composition_conservation]: delta=%.3e\n", delta2);

}


/* ================================================= */

/* compute melt temperature Tm(P) */
static void Tm(struct All_variables *E) {
  int i, n;
  double radius, P, T0, A;

  const int m = 1;
  const double k = 2 * 3.14 * 6.67e-11 / 3
      * 3400 * 3400 * 1e6 * 1e-9 * pow(E->data.radius_km, 2);

  for (i = 1; i <= E->lmesh.nno; i++)
    for (n = 1; n <= E->mush.consts.ncomp; n++) {
      radius = E->sx[m][3][i];
      P = k * (1 - pow(radius, 2));
      T0 = E->mush.consts.melting_point[n];
      A = E->mush.consts.melt_pressure_coff[n];
      E->mush.comp.Tm[n][i] = T0 + A * P;
    }
}

/* partition coefficients K(T, P) */
static double K(struct All_variables *E, double T, int i, int n) {
  double Tm1, K;

  const double T0 = E->mush.consts.melting_point[n];
  const double dS = E->mush.consts.entropy_contrast;
  const double r = E->mush.consts.specific_gas_const[n];

  Tm1 = E->mush.comp.Tm[n][i];
  K = exp(T0 * dS / r * (1 / T - 1 / Tm1));

  return K;
}

/* Newton residual */
static double residual(struct All_variables *E, double T, double phi, int i) {
  int n;
  double Kn, r;

  r = 0;
  for (n = 1; n <= E->mush.consts.ncomp; n++) {
    Kn = K(E, T, i, n);
    r += (1 - Kn) * E->mush.comp.c0[n][i] / (phi + (1 - phi) * Kn);
  }

  return r;
}

/* Newton iteration for T */
static double Newton_iteration_T(struct All_variables *E, double phi, int i) {
  int n, its, flag;
  double T, r, rp, rm, drdT;

  const double r_tol = 1e-8; // tolerance for Newton residual
  const double its_tol = 1000; // maximum number of iterations
  const double eps_T = 5; // temperature perturbation for finite differencing, degrees

  // set starting guess for Tsol
  T = 0;
  for (n = 1; n <= E->mush.consts.ncomp; n++)
    T += E->mush.comp.c0[n][i] * E->mush.comp.Tm[n][i];

  //  get residual for sum(c0^i/K^i) = 1
  r = residual(E, T, phi, i);
  its = 0;
  flag = 1;

  // Newton iteration
  while ((fabs(r) > r_tol) && flag) {
    rp = residual(E, T + eps_T, phi, i);
    rm = residual(E, T - eps_T, phi, i);
    drdT = (rp - rm) / eps_T / 2;
    T -= r / drdT;
    r = residual(E, T, phi, i);

    its += 1;
    if ((its == its_tol) || isinf(T) || isnan(T)) {
      fprintf(stderr, "!!! Newton solver for phi = %.1f, c0=%.4g T=%.4g has not converged !!!\n",
              phi, E->mush.comp.c0[2][i], T);
      T = phi * E->mush.comp.Tliq[i] + (1 - phi) * E->mush.comp.Tsol[i];
      flag = 0;
    }
  }

  return T;
}

/* Newton iteration for phi */
static double Newton_iteration_phi(struct All_variables *E, double T, int i) {
  int n, its;
  double Kn, phi, r, drdphi;

  const double r_tol = 1e-10; // tolerance for Newton residual
  const double its_tol = 100; // maximum number of iterations

  // set starting guess for Tsol
  // phi = (T - E->mush.comp.Tsol[i]) / (E->mush.comp.Tliq[i] - E->mush.comp.Tsol[i]);
  phi = 0;

  //  get residual for sum(c0^i/K^i) = 1
  r = residual(E, T, phi, i);
  its = 0;

  // Newton iteration
  while ((fabs(r) > r_tol) && (its < its_tol)) {
    drdphi = 0;
    for (n = 1; n <= E->mush.consts.ncomp; n++) {
      Kn = K(E, T, i, n);
      drdphi -= E->mush.comp.c0[n][i] * pow((Kn - 1) / (phi + (1 - phi) * Kn), 2);
    }

    phi -= r / drdphi;
    r = residual(E, T, phi, i);

    its += 1;
    if (its == its_tol)
      fprintf(stderr, "!!! Newton solver for equilibrium has not converged !!!\n");
  }

  return phi;
}

/* solid flow */
static void solid_flow(struct All_variables *E, int el, double dc0[9], int n) {
  int i, j;
  float VV[4][9];
  double rtf[4][9];
  double dc0dx[9], dc0dy[9], dc0dz[9], sint[9];
  double vx[9], vy[9], vz[9];
  double c0, Na;
  struct Shape_function PG;

  const int m = 1;
  const int sphere_key = 1;
  const int lev = E->mesh.levmax;
  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int vpts = vpoints[dims];
  const double diff = 0;  // no diffusion term

  velo_from_element(E, VV, m, el, sphere_key);
  get_rtf_at_vpts(E, m, lev, el, rtf);
  pg_shape_fn(E, el, &PG, &(E->gNX[m][el]), VV, rtf, diff, m);

  for (i = 1; i <= vpts; i++) {
    vx[i] = dc0dx[i] = 0.0;
    vy[i] = dc0dy[i] = 0.0;
    vz[i] = dc0dz[i] = 0.0;
    sint[i] = rtf[3][i] / sin(rtf[1][i]);
  }

  for (j = 1; j <= ends; j++) {
    c0 = E->mush.comp.c0[n][E->ien[m][el].node[j]];
    for (i = 1; i <= vpts; i++) {
      dc0dx[i] += c0 * E->gNX[m][el].vpt[GNVXINDEX(0, j, i)] * rtf[3][i];
      dc0dy[i] += c0 * E->gNX[m][el].vpt[GNVXINDEX(1, j, i)] * sint[i];
      dc0dz[i] += c0 * E->gNX[m][el].vpt[GNVXINDEX(2, j, i)];
      Na = E->N.vpt[GNVINDEX(j, i)];
      vx[i] += VV[1][j] * Na;
      vy[i] += VV[2][j] * Na;
      vz[i] += VV[3][j] * Na;
    }
  }

  for (j = 1; j <= ends; j++) {
    dc0[j] = 0;
    for (i = 1; i <= vpts; i++)
      dc0[j] += (vx[i] * dc0dx[i] + vy[i] * dc0dy[i] + vz[i] * dc0dz[i]) *
          PG.vpt[GNVINDEX(j, i)] * E->gDA[m][el].vpt[i];
  }
}