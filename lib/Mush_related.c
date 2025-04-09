#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>
#include <stdio.h>

static void melt_extraction(struct All_variables *E);
static void update_segregation_flux(struct All_variables *E);
static void update_melt_energy(struct All_variables *E);
static void boundary_melt_flux(struct All_variables *E);

void mush_input(struct All_variables *E) {
  int n;
  double z;

  const int m = E->parallel.me;
  const double Myr = 365.25 * 24 * 3600 * 1e6;
  const double phi0 = 0.2;

  /* main */
  input_boolean("mush_version", &(E->mush.mush_version), "off", m);

  /* constants */
  input_double("latent_heat", &(E->mush.consts.latent_heat), "4e5", m);
  input_double("density_contrast", &(E->mush.consts.density_contrast), "500", m);
  input_float("melt_viscosity", &(E->data.melt_viscosity), "1", m);
  input_double("visc_melt_coff", &(E->mush.consts.visc_melt_coff), "26", m);
  input_float("permeability", &(E->data.permeability), "3e-10", m);
  input_double("tensile_strength", &(E->mush.consts.tensile_strength), "5e6", m);
  input_double("extraction_coff", &(E->mush.consts.extraction_coff), "1.4e-5", m);
  input_double("extraction_depth", &(E->mush.consts.extraction_depth), "200", m);
  input_string("background_profile_file", E->control.background_profile_file, "", m);

  /* composition */
  input_int("num_composition", &(n), "1", m);
  E->mush.consts.ncomp = n;
  input_boolean("constant_density", &(E->mush.comp.constant_density), "on", m);
  input_double_vector("solid_density", n, (E->mush.consts.solid_density), m);
  input_double_vector("melt_density", n, (E->mush.consts.melt_density), m);
  input_double_vector("bulk_composition", n, E->mush.consts.bulk_comp, m);
  input_double_vector("melting_point", n, (E->mush.consts.melting_point), m);
  input_double_vector("specific_gas_const", n, (E->mush.consts.specific_gas_const), m);
  input_double_vector("melt_pressure_coff", n, (E->mush.consts.melt_pressure_coff), m);
  input_double("entropy_contrast", &(E->mush.consts.entropy_contrast), "320", m);

  for (n = E->mush.consts.ncomp; n >= 1; n--) {
    E->mush.consts.solid_density[n] = E->mush.consts.solid_density[n - 1];
    E->mush.consts.melt_density[n] = E->mush.consts.melt_density[n - 1];
    E->mush.consts.bulk_comp[n] = E->mush.consts.bulk_comp[n - 1];
    E->mush.consts.melting_point[n] = E->mush.consts.melting_point[n - 1];
    E->mush.consts.specific_gas_const[n] = E->mush.consts.specific_gas_const[n - 1];
    E->mush.consts.melt_pressure_coff[n] = E->mush.consts.melt_pressure_coff[n - 1];
  }

  /* dynamics */
  input_double("convection_developed_timestep",
               &(E->mush.dynamics.convection_developed_timestep), "5e-9", m);
  input_boolean("darcy_migration", &(E->mush.dynamics.darcy_migration), "off", m);
  input_boolean("melt_extraction", &(E->mush.dynamics.dike_extraction), "off", m);
  input_boolean("surface_downwelling", &(E->mush.dynamics.surface_downwelling), "off", m);
  input_boolean("melt_buoyancy", &(E->mush.dynamics.melt_buoyancy), "off", m);
  input_float("extraction_limit", &(E->mush.dynamics.extraction_limit), "0.5", m);
  input_float("downwelling_coff", &(E->mush.dynamics.downwelling_coff), "1", m);
  input_boolean("solve_CMB_temp", &(E->mush.dynamics.solve_CMB_temp), "off", m);
  input_boolean("solve_CMB_comp", &(E->mush.dynamics.solve_CMB_comp), "off", m);
  input_boolean("inhomo_radio", &(E->mush.dynamics.inhomo_radio), "off", m);
  input_double("cpc", &(E->mush.dynamics.Cpc), "1000", m);
  input_boolean("solid_compaction", &(E->mush.dynamics.solid_compaction), "off", m);
  input_boolean("segregation_energy", &(E->mush.dynamics.segregation_energy), "on", m);

  // nondim
  E->mush.consts.St0 = E->mush.consts.latent_heat / (E->data.Cp * E->data.ref_temperature);
  E->mush.consts.B = E->mush.consts.density_contrast /
      (E->data.therm_exp * E->data.res_density * E->data.ref_temperature);
  E->mush.consts.Kq = E->data.permeability * E->mush.consts.density_contrast * E->data.grav_acc *
      E->data.radius_km * 1e3 / (E->data.melt_viscosity * E->data.therm_diff);
  E->mush.consts.KE = Myr / (E->mush.consts.extraction_coff * E->data.ref_viscosity);
  E->mush.consts.KE0 = E->mush.consts.tensile_strength * pow(E->data.radius_km * 1e3, 2) /
      (E->data.ref_viscosity * E->data.therm_diff);
  z = E->data.radius_km / E->mush.consts.extraction_depth;
  E->mush.consts.KE1 = 3 * pow(z, 2) * (pow(phi0, 2) - 0.667 * pow(phi0, 3))
      * E->mush.consts.Kq - 2 * z * E->mush.consts.KE0;

  if (E->parallel.me == 0)
    fprintf(stderr, "(mush nondim) St=%.3g, B=%.3g, Kq=%.3g, KE=%.3g, KE0=%.3g, KE1=%.3g\n",
            E->mush.consts.St0, E->mush.consts.B, E->mush.consts.Kq,
            E->mush.consts.KE, E->mush.consts.KE0, E->mush.consts.KE1);
}

void mush_allocate(struct All_variables *E) {
  int nno = E->lmesh.nno;
  int i, n;

  /* composition */
  E->mush.comp.Tsol = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.comp.Tliq = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.comp.phi = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.comp.drho = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.comp.drho2 = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.comp.St = (double *) malloc((nno + 1) * sizeof(double));
  for (n = 1; n <= E->mush.consts.ncomp; n++) {
    E->mush.comp.Tm[n] = (double *) malloc((nno + 1) * sizeof(double));
    E->mush.comp.c0[n] = (double *) malloc((nno + 1) * sizeof(double));
    E->mush.comp.cl[n] = (double *) malloc((nno + 1) * sizeof(double));
    E->mush.comp.cs[n] = (double *) malloc((nno + 1) * sizeof(double));
  }

  /* dynamics */
  E->mush.dynamics.q = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.dynamics.E = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.dynamics.M = (double *) malloc((nno + 1) * sizeof(double));
  E->mush.dynamics.H = (double *) malloc((nno + 1) * sizeof(double));
  for (n = 1; n <= E->mush.consts.ncomp; n++)
    E->mush.dynamics.Mn[n] = (double *) malloc((nno + 1) * sizeof(double));

  /* init */
  for (i = 1; i <= nno; i++) {
    E->mush.dynamics.E[i] = E->mush.dynamics.M[i] = 0;
    E->mush.dynamics.q[i] = E->mush.dynamics.H[i] = 0;

    E->mush.comp.phi[i] = E->mush.comp.St[i] = 0;
    for (n = 1; n <= E->mush.consts.ncomp; n++) {
      E->mush.comp.c0[n][i] = E->mush.consts.bulk_comp[n];
      E->mush.dynamics.Mn[n][i] = 0;
    }

    // CMB boundary
//    if ((E->parallel.me_loc[3] == 0) && ((i - 1) % E->lmesh.noz == 0)) {
//      E->mush.comp.c0[1][i] = 1;
//      E->mush.comp.c0[2][i] = 0;
//    }
  }

  if (E->mush.dynamics.darcy_migration)
    E->mush.dynamics.convection_developed = -10;
  else
    E->mush.dynamics.convection_developed = 1;
}

/* melt generation, darcy migration and dike extraction */
void process_melting(struct All_variables *E) {
  // update bulk composition
  update_bulk_composition(E);

  // thermodynamic calculation
  thermodynamic_calculation(E);

  // darcy migration
  if (E->mush.dynamics.darcy_migration)
    update_segregation_flux(E);

  // dike extraction
  if (E->mush.dynamics.dike_extraction)
    melt_extraction(E);

  // melt-related energy
  update_melt_energy(E);

  // boundary melt-related heatflux
  boundary_melt_flux(E);
}

/* deal with surface downwelling */
void surface_downwelling(struct All_variables *E) {
  int i, j, isurf;
  double dr;
  float wD;

  const int m = 1;
  const int nno = E->lmesh.nno;
  const int noz = E->lmesh.noz;
  const float k = E->mush.dynamics.downwelling_coff;
  const double late_timescale = 4e-3; // 400 Myr

  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    for (i = 1; i <= nno; i++)
      if (E->mush.dynamics.E[i]) {
        isurf = ((i - 1) / noz + 1) * noz;
        dr = E->sx[m][3][i + 1] - E->sx[m][3][i];
        wD = (float) (k * pow(E->sx[m][3][i], 2) * E->mush.dynamics.E[i] * dr);

        // late state. only local volanism. don't consider downwelling
        if (E->monitor.elapsed_time > late_timescale)
          wD = 0;

        for (j = i; j <= isurf; j++)
          E->sphere.cap[m].V[3][j] -= wD / (float) pow(E->sx[m][3][j], 2);
      }
}


/* ================================================== */

/* segregation flux q */
static void update_segregation_flux(struct All_variables *E) {
  int i;
  double phi, drho;
  const double Kq = E->mush.consts.Kq;
//  const double phic = 0.4;
//  const double qmax = Kq * pow(phic, 3);

  if (E->mush.dynamics.convection_developed > 0)
    for (i = 1; i <= E->lmesh.nno; i++) {
      phi = E->mush.comp.phi[i];
      // drho = E->mush.comp.drho[i];
      drho = 1;
      // E->mush.dynamics.q[i] = Kq * drho * pow(phi, 3) * pow((1 - phi), 2);
      // E->mush.dynamics.q[i] = min(E->mush.dynamics.q[i], qmax);
      E->mush.dynamics.q[i] = Kq * drho * pow(phi, 3);
    }
  else if (E->advection.timestep < E->mush.dynamics.convection_developed_timestep)
    E->mush.dynamics.convection_developed++;
}

/* melt extraction (only for the top caps) */
static void melt_extraction(struct All_variables *E) {
  int i, j, el, isurf, n, node1, node2;
  double *E1[NCS], e1, dqdr, phi, z, zeta, Pc, dr, dphi, mhat, rhol, rhos;

  const int m = 1;
  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int noz = E->lmesh.noz;
  const int nc = E->mush.consts.ncomp;
  const int lev = E->mesh.levmax;
  const float dt = E->advection.timestep;
  const double KE0 = E->mush.consts.KE0;
  const double KE = E->mush.consts.KE;
  const double KE1 = E->mush.consts.KE1;
  const double limit = E->mush.dynamics.extraction_limit;
  const double visE = E->viscosity.E[0];
  const double visT = E->viscosity.T[0];
  const double rho0 = E->data.density;
  const double drho0 = E->mush.consts.density_contrast;
  const double phic = 0.4;
  const double weight = 0.25;

  E1[m] = (double *) malloc((nno + 1) * sizeof(double));
  for (i = 1; i <= nno; i++) {
    E1[m][i] = E->mush.dynamics.M[i] = 0;
    for (n = 1; n <= nc; n++)
      E->mush.dynamics.Mn[n][i] = 0;
  }

  // melt extraction rate E
  for (el = 1; el <= nel; el++) {
    dr = E->eco[m][el].size[3];
    for (j = 5; j <= 8; j++) {
      node2 = E->ien[m][el].node[j];
      node1 = E->ien[m][el].node[j - 4];
      phi = E->mush.comp.phi[node1];

      // melt density should be smaller than crust density
      rhol = rho0 - drho0 *
          (E->mush.comp.drho2[node1] + (1 - phi) * E->mush.comp.drho[node1]);
      rhos = rho0 - drho0 * E->mush.comp.drho2[((node1 - 1) / noz + 1) * noz];

      rhol = rhos - 0.1; // ignore it now

      if (phi && rhol < rhos) {
        z = 1 - E->sx[m][3][node1];
        dqdr = (E->mush.dynamics.q[node2] - E->mush.dynamics.q[node1]) / dr;
        dphi = fabs(E->mush.comp.phi[node1] - E->mush.comp.phi[node2]);

        // late state. only local melting. make case10 converge.
        if (E->monitor.elapsed_time > 45e-3)
          phi = min(phi, 0.02);

        zeta = (1 - min(phi, phic)) * exp(visE * (1 / E->T[m][node1] - 1 / visT)) / phi;
        Pc = -zeta * dqdr;
        e1 = max(0, (Pc - KE1 * z - KE0) / (KE + zeta));
        e1 = min(e1, dphi / dt * limit);
        E1[m][node1] += e1 * weight;

        // output debug
//        if (((E->parallel.me == 1) && (el <= 32) && (el >= 1) && (j == 5)))
//          fprintf(stderr, "dike[%d]: zeta=%6.4g dqdr=%8.2g Pc=%8.2g Bz=%8.2g\n",
//                  el, zeta, dqdr, Pc, B * z);
      }
    }
  }

  (E->exchange_node_d)(E, E1, lev);
  for (i = 1; i <= nno; i++)
    E->mush.dynamics.E[i] = E1[m][i];

  // dike emplacement M
  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
    for (i = 1; i <= nno; i++)
      if (E1[m][i]) {
        isurf = ((i - 1) / noz + 1) * noz;
        // mhat = E1[m][i] / (1 - E->sx[m][3][i]);
        if (i % E->lmesh.noz == 1)
          dr = (E->sx[m][3][i + 1] - E->sx[m][3][i]);
        else
          dr = (E->sx[m][3][i] - E->sx[m][3][i - 1]);
//        for (j = i + 1; j <= isurf; j++) {
//          E->mush.dynamics.M[j] += mhat * dr;
//          for (n = 1; n <= nc; n++)
//            E->mush.dynamics.Mn[n][j] += mhat * E->mush.comp.cl[n][i] * dr;
//        }
        E->mush.dynamics.M[isurf] += pow(E->sx[m][3][i], 2) * E1[m][i] * dr;
        for (n = 1; n <= nc; n++)
          E->mush.dynamics.Mn[n][isurf] +=
              pow(E->sx[m][3][i], 2) * E1[m][i] * (E->mush.comp.cl[n][i] - E->mush.comp.c0[n][i]) * dr;
      }
  }

  free(E1[m]);
}

/* update melt energy */
static void update_melt_energy(struct All_variables *E) {
  int el, i, j, node;
  double phi, q, dTdr, dqdr;
  double Na, dNadr, dV;
  double heatflux, massflux;
  double *H[NCS];

  const int m = 1;
  const int nno = E->lmesh.nno;
  const int lev = E->mesh.levmax;
  const double L = E->mush.consts.St0;

  H[m] = (double *) malloc((nno + 1) * sizeof(double));
  for (i = 1; i <= nno; i++)
    H[m][i] = 0;

  for (el = 1; el <= E->lmesh.nel; el++) {
    heatflux = massflux = 0;
    for (i = 1; i <= 8; i++) {
      phi = q = dqdr = dTdr = 0;
      for (j = 1; j <= 8; j++) {
        Na = E->N.vpt[GNVINDEX(j, i)];
        dNadr = E->gNX[m][el].vpt[GNVXINDEX(2, j, i)];
        node = E->ien[m][el].node[j];
        phi += E->mush.comp.phi[node] * Na;
        q += E->mush.dynamics.q[node] * Na;
        dTdr += E->T[m][node] * dNadr;
        dqdr += E->mush.dynamics.q[node] * dNadr;
      }
      dV = E->gDA[m][el].vpt[i] / E->eco[m][el].area;
      heatflux += q * dTdr * dV;
      massflux += L * (1 - phi) * dqdr * dV;
    }
    for (j = 5; j <= 8; j++)
      H[m][E->ien[m][el].node[j]] += (heatflux + massflux) * 0.25;
  }
  (E->exchange_node_d)(E, H, lev);
  for (i = 1; i <= nno; i++) {
    if (E->mush.dynamics.segregation_energy)
      E->mush.dynamics.H[i] = -H[m][i] - L * (1 - E->mush.comp.phi[i]) * E->mush.dynamics.E[i];
    else
      E->mush.dynamics.H[i] = -L * (1 - E->mush.comp.phi[i]) * E->mush.dynamics.E[i];
  }

  free(H[m]);
}

/* melt flux (ZYZ-220412) */
static void boundary_melt_flux(struct All_variables *E) {
  int e, el, i, j, node;
  double q, T, qT;
  double Na, dV;
  float *sum_h;
  float *flux[NCS];

  const double L = E->mush.consts.St0;
  const int vpts = vpoints[E->mesh.nsd];
  const int ends = enodes[E->mesh.nsd];
  const int nno = E->lmesh.nno;
  const int lev = E->mesh.levmax;
  const int m = 1;

  sum_h = (float *) malloc((5) * sizeof(float));
  for (i = 0; i <= 4; i++)
    sum_h[i] = 0;

  flux[m] = (float *) malloc((1 + nno) * sizeof(float));
  for (i = 1; i <= nno; i++)
    flux[m][i] = 0;

  for (e = 1; e <= E->lmesh.nel; e++) {
    qT = 0;
    for (i = 1; i <= vpts; i++) {
      q = T = 0;
      for (j = 1; j <= ends; j++) {
        node = E->ien[m][e].node[j];
        Na = E->N.vpt[GNVINDEX(j, i)];
        q += E->mush.dynamics.q[node] * Na;
        T += E->T[m][node] * Na;
      }
      dV = E->gDA[m][e].vpt[i] / E->eco[m][e].area;
      qT += q * (T + L) * dV;
    }
    for (j = 1; j <= ends; j++)
      flux[m][E->ien[m][e].node[j]] += (float) qT * E->TWW[lev][m][e].node[j];
  }

  (E->exchange_node_f)(E, flux, lev);

  for (i = 1; i <= nno; i++)
    flux[m][i] *= (float) (E->MASS[lev][m][i] / (1 + E->mush.comp.St[i]));

  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1)
    for (i = 1; i <= E->lmesh.nsf; i++)
      E->slice.smflux[m][i] = flux[m][E->surf_node[m][i]] +
          (float) ((E->T[m][E->surf_node[m][i]] + L) * E->mush.dynamics.M[E->surf_node[m][i]]);

  if (E->parallel.me_loc[3] == 0)
    for (i = 1; i <= E->lmesh.nsf; i++)
      E->slice.bmflux[m][i] = flux[m][E->surf_node[m][i] - E->lmesh.noz + 1];

  for (e = 1; e <= E->lmesh.snel; e++) {
    qT = (float) ((E->slice.smflux[m][E->sien[m][e].node[1]] +
        E->slice.smflux[m][E->sien[m][e].node[2]] +
        E->slice.smflux[m][E->sien[m][e].node[3]] +
        E->slice.smflux[m][E->sien[m][e].node[4]]) *
        0.25);
    el = e * E->lmesh.elz;
    sum_h[0] += (float) qT * E->eco[m][el].area;
    sum_h[1] += E->eco[m][el].area;

    qT = (float) ((E->slice.bmflux[m][E->sien[m][e].node[1]] +
        E->slice.bmflux[m][E->sien[m][e].node[2]] +
        E->slice.bmflux[m][E->sien[m][e].node[3]] +
        E->slice.bmflux[m][E->sien[m][e].node[4]]) *
        0.25);
    el = (e - 1) * E->lmesh.elz + 1;
    sum_h[2] += (float) qT * E->eco[m][el].area;
    sum_h[3] += E->eco[m][el].area;
  }

  sum_across_surface(E, sum_h, 4);

  if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
    sum_h[0] = sum_h[0] / sum_h[1];
    E->sphere.cap[m].melt_flux = sum_h[0];
//    if (E->parallel.me == E->parallel.nprocz - 1)
//      fprintf(stderr, "surface melt flux= %.4g\n", sum_h[0]);
  }

  if (E->parallel.me_loc[3] == 0) {
    sum_h[2] = sum_h[2] / sum_h[3];
    E->sphere.cap[m].melt_flux = sum_h[2];
//    if (E->parallel.me == 0)
//      fprintf(stderr, "bottom melt flux= %.4g\n", sum_h[2]);
  }

  free((void *) flux[m]);
  free((void *) sum_h);
}