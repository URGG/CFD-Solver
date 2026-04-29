#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_splinalg.h>
#include "raylib.h"

#define NX 40
#define NY 40

#define IMAX (NX + 2)
#define JMAX (NY + 2)

typedef struct {
    double p[IMAX * JMAX];
    double u[IMAX * JMAX];
    double v[IMAX * JMAX];
    double us[IMAX * JMAX];
    double vs[IMAX * JMAX];
} staggered_grid_t;

static size_t IX(size_t i, size_t j) {
    return j * IMAX + i;
}

static size_t GSL_IX(size_t i, size_t j) {
    return (j - 1) * NX + (i - 1);
}

static void apply_boundary_conditions(staggered_grid_t *grid) {
    for (size_t j = 0; j < JMAX; ++j) {
        grid->u[IX(1, j)] = 0.0;
        grid->u[IX(NX + 1, j)] = 0.0;
        grid->v[IX(0, j)] = -grid->v[IX(1, j)];
        grid->v[IX(NX + 1, j)] = -grid->v[IX(NX, j)];

        grid->us[IX(1, j)] = 0.0;
        grid->us[IX(NX + 1, j)] = 0.0;
    }
    for (size_t i = 0; i < IMAX; ++i) {
        grid->v[IX(i, 1)] = 0.0;
        grid->v[IX(i, NY + 1)] = 0.0;
        grid->u[IX(i, 0)] = -grid->u[IX(i, 1)];
        grid->u[IX(i, NY + 1)] = 2.0 * 1.0 - grid->u[IX(i, NY)];

        grid->vs[IX(i, 1)] = 0.0;
        grid->vs[IX(i, NY + 1)] = 0.0;
    }
}

static void compute_tentative_velocity(staggered_grid_t *grid, double dt, double nu, double dx, double dy) {
    const double dxi = 1.0 / dx;
    const double dyi = 1.0 / dy;
    const double dxi2 = dxi * dxi;
    const double dyi2 = dyi * dyi;

    for (size_t j = 1; j <= NY; ++j) {
        for (size_t i = 2; i <= NX; ++i) {
            double v_here = 0.25 * (grid->v[IX(i-1, j)] + grid->v[IX(i-1, j+1)] + grid->v[IX(i, j)] + grid->v[IX(i, j+1)]);
            double u_center = grid->u[IX(i, j)];

            double diff_x = nu * (grid->u[IX(i-1, j)] - 2.0 * u_center + grid->u[IX(i+1, j)]) * dxi2;
            double diff_y = nu * (grid->u[IX(i, j-1)] - 2.0 * u_center + grid->u[IX(i, j+1)]) * dyi2;

            double adv_x = u_center * (grid->u[IX(i+1, j)] - grid->u[IX(i-1, j)]) * 0.5 * dxi;
            double adv_y = v_here * (grid->u[IX(i, j+1)] - grid->u[IX(i, j-1)]) * 0.5 * dyi;

            grid->us[IX(i, j)] = u_center + dt * (diff_x + diff_y - adv_x - adv_y);
        }
    }

    for (size_t j = 2; j <= NY; ++j) {
        for (size_t i = 1; i <= NX; ++i) {
            double u_here = 0.25 * (grid->u[IX(i, j-1)] + grid->u[IX(i, j)] + grid->u[IX(i+1, j-1)] + grid->u[IX(i+1, j)]);
            double v_center = grid->v[IX(i, j)];

            double diff_x = nu * (grid->v[IX(i-1, j)] - 2.0 * v_center + grid->v[IX(i+1, j)]) * dxi2;
            double diff_y = nu * (grid->v[IX(i, j-1)] - 2.0 * v_center + grid->v[IX(i, j+1)]) * dyi2;

            double adv_x = u_here * (grid->v[IX(i+1, j)] - grid->v[IX(i-1, j)]) * 0.5 * dxi;
            double adv_y = v_center * (grid->v[IX(i, j+1)] - grid->v[IX(i, j-1)]) * 0.5 * dyi;

            grid->vs[IX(i, j)] = v_center + dt * (diff_x + diff_y - adv_x - adv_y);
        }
    }
}

static void build_laplacian(gsl_spmatrix *A, double dx, double dy) {
    double dxi2 = 1.0 / (dx * dx);
    double dyi2 = 1.0 / (dy * dy);

    for (size_t j = 1; j <= NY; ++j) {
        for (size_t i = 1; i <= NX; ++i) {
            size_t row = GSL_IX(i, j);
            double diag = 0.0;

            if (i > 1) { gsl_spmatrix_set(A, row, GSL_IX(i - 1, j), -dxi2); diag += dxi2; }
            if (i < NX) { gsl_spmatrix_set(A, row, GSL_IX(i + 1, j), -dxi2); diag += dxi2; }
            if (j > 1) { gsl_spmatrix_set(A, row, GSL_IX(i, j - 1), -dyi2); diag += dyi2; }
            if (j < NY) { gsl_spmatrix_set(A, row, GSL_IX(i, j + 1), -dyi2); diag += dyi2; }

            if (i == 1 && j == 1) {
                diag = 1.0;
                gsl_spmatrix_set(A, row, GSL_IX(i + 1, j), 0.0);
                gsl_spmatrix_set(A, row, GSL_IX(i, j + 1), 0.0);
            }

            gsl_spmatrix_set(A, row, row, diag);
        }
    }
}

static void solve_pressure(staggered_grid_t *grid, gsl_spmatrix *A_csr, gsl_splinalg_itersolve *solver, double dt, double rho, double dx, double dy) {
    gsl_vector *b = gsl_vector_alloc(NX * NY);
    gsl_vector *x = gsl_vector_alloc(NX * NY);

    for (size_t j = 1; j <= NY; ++j) {
        for (size_t i = 1; i <= NX; ++i) {
            size_t row = GSL_IX(i, j);
            double div = ((grid->us[IX(i + 1, j)] - grid->us[IX(i, j)]) / dx) +
                         ((grid->vs[IX(i, j + 1)] - grid->vs[IX(i, j)]) / dy);

            gsl_vector_set(b, row, -(rho / dt) * div);
            gsl_vector_set(x, row, grid->p[IX(i, j)]);
        }
    }

    gsl_vector_set(b, GSL_IX(1, 1), 0.0);

    int status;
    do {
        status = gsl_splinalg_itersolve_iterate(A_csr, b, 1e-5, x, solver);
    } while (status == GSL_CONTINUE);

    for (size_t j = 1; j <= NY; ++j) {
        for (size_t i = 1; i <= NX; ++i) {
            grid->p[IX(i, j)] = gsl_vector_get(x, GSL_IX(i, j));
        }
    }

    gsl_vector_free(b);
    gsl_vector_free(x);
}

static void apply_corrector(staggered_grid_t *grid, double dt, double rho, double dx, double dy) {
    double dxi = 1.0 / dx;
    double dyi = 1.0 / dy;

    for (size_t j = 1; j <= NY; ++j) {
        for (size_t i = 2; i <= NX; ++i) {
            grid->u[IX(i, j)] = grid->us[IX(i, j)] - (dt / rho) * (grid->p[IX(i, j)] - grid->p[IX(i-1, j)]) * dxi;
        }
    }
    for (size_t j = 2; j <= NY; ++j) {
        for (size_t i = 1; i <= NX; ++i) {
            grid->v[IX(i, j)] = grid->vs[IX(i, j)] - (dt / rho) * (grid->p[IX(i, j)] - grid->p[IX(i, j-1)]) * dyi;
        }
    }
}

int main(void) {
    staggered_grid_t *grid = calloc(1, sizeof(staggered_grid_t));

    double dt = 0.001;
    double rho = 1.0;
    double nu = 0.01;
    double dx = 1.0 / NX;
    double dy = 1.0 / NY;

    gsl_spmatrix *A = gsl_spmatrix_alloc(NX * NY, NX * NY);
    build_laplacian(A, dx, dy);
    gsl_spmatrix *A_csr = gsl_spmatrix_crs(A);
    gsl_splinalg_itersolve *solver = gsl_splinalg_itersolve_alloc(gsl_splinalg_itersolve_gmres, NX * NY, 0);

    const int cellSize = 20;
    InitWindow(NX * cellSize, NY * cellSize, "Navier-Stokes Lid-Driven Cavity");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        for (int k = 0; k < 5; k++) {
            apply_boundary_conditions(grid);
            compute_tentative_velocity(grid, dt, nu, dx, dy);
            apply_boundary_conditions(grid);
            solve_pressure(grid, A_csr, solver, dt, rho, dx, dy);
            apply_corrector(grid, dt, rho, dx, dy);
            apply_boundary_conditions(grid);
        }

        BeginDrawing();
        ClearBackground(BLACK);

        for (size_t j = 1; j <= NY; ++j) {
            for (size_t i = 1; i <= NX; ++i) {
                double u_c = 0.5 * (grid->u[IX(i, j)] + grid->u[IX(i+1, j)]);
                double v_c = 0.5 * (grid->v[IX(i, j)] + grid->v[IX(i, j+1)]);

                double speed = sqrt(u_c * u_c + v_c * v_c);
                if (speed > 1.0) speed = 1.0;

                unsigned char c = (unsigned char)(speed * 255.0);
                DrawRectangle((i - 1) * cellSize, (NY - j) * cellSize, cellSize, cellSize, (Color){c, c, 255, 255});
            }
        }

        DrawFPS(10, 10);
        EndDrawing();
    }

    gsl_splinalg_itersolve_free(solver);
    gsl_spmatrix_free(A_csr);
    gsl_spmatrix_free(A);
    free(grid);
    CloseWindow();

    return EXIT_SUCCESS;
}