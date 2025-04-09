/* Functions relating to the building and use of mesh locations ... */


#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>

static void compute_sphereh_table(struct All_variables *);

/*   ======================================================================
     ======================================================================  */

void set_sphere_harmonics(E)
    struct All_variables *E;

{
  int m, node, ll, mm, i, j;

  i = 0;
  for (ll = 0; ll <= E->output.llmax; ll++)
    for (mm = 0; mm <= ll; mm++) {
      E->sphere.hindex[ll][mm] = i;
      i++;
    }

  E->sphere.hindice = i;

  /* spherical harmonic coeff (0=cos, 1=sin)
     for surface topo, cmb topo and geoid */
  for (i = 0; i <= 1; i++) {
    E->sphere.harm_geoid[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    E->sphere.harm_geoid_from_bncy[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    E->sphere.harm_geoid_from_bncy_botm[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    E->sphere.harm_geoid_from_tpgt[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    E->sphere.harm_geoid_from_tpgb[i] = (float *) malloc(E->sphere.hindice * sizeof(float));

    E->sphere.harm_tpgt[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    E->sphere.harm_tpgb[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
  }

  compute_sphereh_table(E);

  return;
}

/* =====================================
   Generalized Legendre polynomials
   =====================================*/
double modified_plgndr_a(int l, int m, double t) {
  int i, ll;
  double x, fact1, fact2, fact, pll, pmm, pmmp1, somx2, plgndr;
  const double three = 3.0;
  const double two = 2.0;
  const double one = 1.0;

  x = cos(t);
  pmm = one;
  if (m > 0) {
    somx2 = sqrt((one - x) * (one + x));
    fact1 = three;
    fact2 = two;
    for (i = 1; i <= m; i++) {
      fact = sqrt(fact1 / fact2);
      pmm = -pmm * fact * somx2;
      fact1 += two;
      fact2 += two;
    }
  }

  if (l == m)
    plgndr = pmm;
  else {
    pmmp1 = x * sqrt(two * m + three) * pmm;
    if (l == m + 1)
      plgndr = pmmp1;
    else {
      for (ll = m + 2; ll <= l; ll++) {
        fact1 = sqrt((4.0 * ll * ll - one) * (double) (ll - m) / (double) (ll + m));
        fact2 = sqrt((2.0 * ll + one) * (ll - m) * (ll + m - one) * (ll - m - one)
                         / (double) ((two * ll - three) * (ll + m)));
        pll = (x * fact1 * pmmp1 - fact2 * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      plgndr = pll;
    }
  }

  plgndr /= sqrt(4.0 * M_PI);

  if (m != 0)
    plgndr *= sqrt(two);

  return plgndr;
}

/* =========================================================
   expand the field TG into spherical harmonics
   ========================================================= */
void sphere_expansion(E, TG, sphc, sphs)
    struct All_variables *E;
    float **TG, *sphc, *sphs;
{
  int el, nint, d, p, i, m, j, es, mm, ll, rand();
  void sum_across_surf_sph1();

  for (i = 0; i < E->sphere.hindice; i++) {
    sphc[i] = 0.0;
    sphs[i] = 0.0;
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++)
    for (es = 1; es <= E->lmesh.snel; es++) {

      for (ll = 0; ll <= E->output.llmax; ll++)
        for (mm = 0; mm <= ll; mm++) {

          p = E->sphere.hindex[ll][mm];

          for (nint = 1; nint <= onedvpoints[E->mesh.nsd]; nint++) {
            for (d = 1; d <= onedvpoints[E->mesh.nsd]; d++) {
              j = E->sien[m][es].node[d];
              sphc[p] += TG[m][E->sien[m][es].node[d]]
                  * E->sphere.tablesplm[m][j][p]
                  * E->sphere.tablescosf[m][j][mm]
                  * E->M.vpt[GMVINDEX(d, nint)]
                  * E->surf_det[m][nint][es];
              sphs[p] += TG[m][E->sien[m][es].node[d]]
                  * E->sphere.tablesplm[m][j][p]
                  * E->sphere.tablessinf[m][j][mm]
                  * E->M.vpt[GMVINDEX(d, nint)]
                  * E->surf_det[m][nint][es];
            }
          }

        }       /* end for ll and mm  */

    }

  sum_across_surf_sph1(E, sphc, sphs);

  return;
}

void debug_sphere_expansion(struct All_variables *E) {
  /* expand temperature field (which should be a sph. harm. load)
   * and output the expansion coeff. to stderr
   */
  int m, i, j, k, p, node;
  int ll, mm;
  float *TT[NCS], *sph_harm[2];

  for (m = 1; m <= E->sphere.caps_per_proc; m++)
    TT[m] = (float *) malloc((E->lmesh.nsf + 1) * sizeof(float));

  /* sin coeff */
  sph_harm[0] = (float *) malloc(E->sphere.hindice * sizeof(float));
  /* cos coeff */
  sph_harm[1] = (float *) malloc(E->sphere.hindice * sizeof(float));

  for (k = 1; k <= E->lmesh.noz; k++) {
    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      for (i = 1; i <= E->lmesh.noy; i++)
        for (j = 1; j <= E->lmesh.nox; j++) {
          node = k + (j - 1) * E->lmesh.noz + (i - 1) * E->lmesh.nox * E->lmesh.noz;
          p = j + (i - 1) * E->lmesh.nox;
          TT[m][p] = E->T[m][node];
        }

    /* expand TT into spherical harmonics */
    sphere_expansion(E, TT, sph_harm[0], sph_harm[1]);

    /* only the first nprocz CPU needs output */
    if (E->parallel.me < E->parallel.nprocz) {
      for (ll = 0; ll <= E->output.llmax; ll++)
        for (mm = 0; mm <= ll; mm++) {
          p = E->sphere.hindex[ll][mm];
          fprintf(stderr, "T expanded layer=%d ll=%d mm=%d -- %12g %12g\n",
                  k + E->lmesh.nzs - 1, ll, mm,
                  sph_harm[0][p], sph_harm[1][p]);
        }
    }
  }

  return;
}


/* ==================================================*/
/* ==================================================*/
static void compute_sphereh_table(E)
    struct All_variables *E;
{
  double modified_plgndr_a();

  int m, node, ll, mm, i, j, p;
  double t, f, mmf;

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    E->sphere.tablesplm[m] = (double **) malloc((E->lmesh.nsf + 1) * sizeof(double *));
    E->sphere.tablescosf[m] = (double **) malloc((E->lmesh.nsf + 1) * sizeof(double *));
    E->sphere.tablessinf[m] = (double **) malloc((E->lmesh.nsf + 1) * sizeof(double *));

    for (i = 1; i <= E->lmesh.nsf; i++) {
      E->sphere.tablesplm[m][i] = (double *) malloc((E->sphere.hindice) * sizeof(double));
      E->sphere.tablescosf[m][i] = (double *) malloc((E->output.llmax + 1) * sizeof(double));
      E->sphere.tablessinf[m][i] = (double *) malloc((E->output.llmax + 1) * sizeof(double));
    }
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    for (j = 1; j <= E->lmesh.nsf; j++) {
      node = j * E->lmesh.noz;
      f = E->sx[m][2][node];
      t = E->sx[m][1][node];
      for (mm = 0; mm <= E->output.llmax; mm++) {
        mmf = (double) (mm) * f;
        E->sphere.tablescosf[m][j][mm] = cos(mmf);
        E->sphere.tablessinf[m][j][mm] = sin(mmf);
      }

      for (ll = 0; ll <= E->output.llmax; ll++)
        for (mm = 0; mm <= ll; mm++) {
          p = E->sphere.hindex[ll][mm];
          E->sphere.tablesplm[m][j][p] = modified_plgndr_a(ll, mm, t);
        }
    }
  }
}

/* output spherical harmonics (ZYZ-220323) */
void output_S_H(struct All_variables *E) {
  int ll, mm, m, i, j, k, snode, node, lev;
  float x, y, *TG[NCS], *TGC[NCS], *power[100], *powerC[100], *sph_harm[2], *sph_harmC[2];
  FILE *fp1;
  char output_file[255];

  lev = E->mesh.levmax;
  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    TG[m] = (float *) malloc((E->lmesh.nsf + 1) * sizeof(float));
    TGC[m] = (float *) malloc((E->lmesh.nsf + 1) * sizeof(float));
  }

  for (i = 0; i <= 1; i++) {
    sph_harm[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
    sph_harmC[i] = (float *) malloc(E->sphere.hindice * sizeof(float));
  }

  for (ll = 0; ll <= E->output.llmax; ll++) {
    power[ll] = (float *) malloc((E->lmesh.noz + 1) * sizeof(float));
    powerC[ll] = (float *) malloc((E->lmesh.noz + 1) * sizeof(float));
  }

  for (i = 1; i <= E->lmesh.noz; i++) {
    for (m = 1; m <= E->sphere.caps_per_proc; m++)
      for (j = 1; j <= E->lmesh.nox; j++)
        for (k = 1; k <= E->lmesh.noy; k++) {
          snode = j + (k - 1) * E->lmesh.nox;
          node = i + (j - 1) * E->lmesh.noz + (k - 1) * E->lmesh.nox * E->lmesh.noz;
          TG[m][snode] = E->T[m][node];
          if (E->mush.mush_version)
            TGC[m][snode] = E->mush.comp.phi[node];
        }

    sphere_expansion(E, TG, sph_harm[0], sph_harm[1]); // temperature
    if (E->mush.mush_version)
      sphere_expansion(E, TGC, sph_harmC[0], sph_harmC[1]); // melt fraction

    for (ll = 0; ll <= E->output.llmax; ll++) {
      x = y = 0;
      for (mm = 0; mm <= ll; mm++) {
        k = E->sphere.hindex[ll][mm];
        x += sph_harm[0][k] * sph_harm[0][k]
            + sph_harm[1][k] * sph_harm[1][k];
        y += sph_harmC[0][k] * sph_harmC[0][k]
            + sph_harmC[1][k] * sph_harmC[1][k];
      }
      power[ll][i] = sqrt(x);
      powerC[ll][i] = sqrt(y);
    }
  }

  if (E->parallel.me < E->parallel.nprocz) {
    sprintf(output_file,
            "%s.power_r.%d.%d",
            E->control.data_file,
            E->parallel.me,
            E->monitor.solution_cycles);
    fp1 = fopen(output_file, "w");

    for (ll = 0; ll <= E->output.llmax; ll++) {
      for (i = 1; i <= E->lmesh.noz; i++) {
        fprintf(fp1, "%d %g %g", ll, E->sphere.R[lev][i], power[ll][i]);
        if (E->mush.mush_version)
          fprintf(fp1, " %g", powerC[ll][i]);
        fprintf(fp1, "\n");
      }
    }
    fclose(fp1);
  }

  for (m = 1; m <= E->sphere.caps_per_proc; m++) {
    free((void *) TG[m]);
    free((void *) TGC[m]);
  }
  for (ll = 0; ll <= E->output.llmax; ll++) {
    free((void *) power[ll]);
    free((void *) powerC[ll]);
  }

  free((void *) sph_harm[0]);
  free((void *) sph_harm[1]);
  free((void *) sph_harmC[0]);
  free((void *) sph_harmC[1]);
}


