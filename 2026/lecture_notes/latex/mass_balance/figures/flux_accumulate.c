/*
 * Forward flowline flux integration for balance velocity.
 *
 * For each ice cell, trace a flowline downstream along the velocity
 * field and deposit that cell's SMB flux at every grid cell the
 * flowline passes through.  Convergence is handled naturally: where
 * N flowlines merge into an ice stream, N× the flux accumulates.
 *
 * Compile:  cc -arch x86_64 -O3 -shared -fPIC -o flux_accumulate.so flux_accumulate.c -lm
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

void flowline_integrate(
    double       *total_flux, /* (nrows*ncols) output — zeroed on entry  */
    const double *smb_flux,   /* (nrows*ncols) ȧ·dx² at each cell       */
    const double *vx,         /* (nrows*ncols) velocity x (+east=+col)   */
    const double *vy,         /* (nrows*ncols) velocity y (+north=-row)  */
    const double *spd,        /* (nrows*ncols) speed                     */
    const int    *valid_mask, /* (nrows*ncols) 1=valid ice, 0=not        */
    const int    *seeds,      /* (n_seeds) flat indices of seed cells    */
    int n_seeds, int nrows, int ncols, int max_steps)
{
    int ncells = nrows * ncols;
    /* Per-cell visit stamp: avoid double-depositing within one flowline */
    int *stamp = (int *)calloc(ncells, sizeof(int));

    for (int s = 0; s < n_seeds; s++) {
        int seed = seeds[s];
        double my_flux = smb_flux[seed];
        if (my_flux <= 0.0) continue;

        int cur_stamp = s + 1;   /* unique stamp for this flowline */

        /* Start at cell centre */
        double pr = (double)(seed / ncols) + 0.5;
        double pc = (double)(seed % ncols) + 0.5;

        for (int step = 0; step < max_steps; step++) {
            int r = (int)pr;
            int c = (int)pc;
            if (r < 0 || r >= nrows || c < 0 || c >= ncols) break;
            int idx = r * ncols + c;
            if (!valid_mask[idx]) break;

            /* Deposit flux (once per cell per flowline) */
            if (stamp[idx] != cur_stamp) {
                stamp[idx] = cur_stamp;
                total_flux[idx] += my_flux;
            }

            /* Bilinear interpolation of unit velocity at (pr, pc) */
            double fr = pr - r;   /* fractional row within cell */
            double fc = pc - c;   /* fractional col within cell */

            /* Clamp neighbour indices */
            int r1 = (r + 1 < nrows) ? r + 1 : r;
            int c1 = (c + 1 < ncols) ? c + 1 : c;

            double ux00 = (spd[r *ncols+c ] > 0.01) ? vx[r *ncols+c ]/spd[r *ncols+c ] : 0.0;
            double ux10 = (spd[r1*ncols+c ] > 0.01) ? vx[r1*ncols+c ]/spd[r1*ncols+c ] : 0.0;
            double ux01 = (spd[r *ncols+c1] > 0.01) ? vx[r *ncols+c1]/spd[r *ncols+c1] : 0.0;
            double ux11 = (spd[r1*ncols+c1] > 0.01) ? vx[r1*ncols+c1]/spd[r1*ncols+c1] : 0.0;

            double uy00 = (spd[r *ncols+c ] > 0.01) ? vy[r *ncols+c ]/spd[r *ncols+c ] : 0.0;
            double uy10 = (spd[r1*ncols+c ] > 0.01) ? vy[r1*ncols+c ]/spd[r1*ncols+c ] : 0.0;
            double uy01 = (spd[r *ncols+c1] > 0.01) ? vy[r *ncols+c1]/spd[r *ncols+c1] : 0.0;
            double uy11 = (spd[r1*ncols+c1] > 0.01) ? vy[r1*ncols+c1]/spd[r1*ncols+c1] : 0.0;

            double ux_interp = ux00*(1-fr)*(1-fc) + ux10*fr*(1-fc)
                             + ux01*(1-fr)*fc     + ux11*fr*fc;
            double uy_interp = uy00*(1-fr)*(1-fc) + uy10*fr*(1-fc)
                             + uy01*(1-fr)*fc     + uy11*fr*fc;

            double mag = sqrt(ux_interp*ux_interp + uy_interp*uy_interp);
            if (mag < 1e-6) break;  /* stagnation — stop */
            ux_interp /= mag;
            uy_interp /= mag;

            /* Advance one cell-width step along the velocity */
            pc += ux_interp;        /* +east  = +col */
            pr -= uy_interp;        /* +north = -row */
        }
    }

    free(stamp);
}
