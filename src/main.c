#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#define NX 41
#define NY 41
#define DX (1.0 / (double)(NX - 1))
#define DY (1.0 / (double)(NY - 1))
#define DIFFUSIVITY 0.05
#define SAFETY 0.20
#define STEPS 800
#define REPORT_EVERY 100

typedef struct scalar_field {
    double values[NX * NY];
} scalar_field_t;

static size_t cell_index(size_t i, size_t j) {
    return j * NX + i;
}

static void field_fill(scalar_field_t *field, double value) {
    for (size_t k = 0; k < NX * NY; ++k) {
        field->values[k] = value;
    }
}

static void apply_boundary_conditions(scalar_field_t *field) {
    for (size_t i = 0; i < NX; ++i) {
        field->values[cell_index(i, 0)] = 0.0;
        field->values[cell_index(i, NY - 1)] = 0.0;
    }

    for (size_t j = 0; j < NY; ++j) {
        field->values[cell_index(0, j)] = 0.0;
        field->values[cell_index(NX - 1, j)] = 0.0;
    }

    for (size_t j = NY / 3; j < (2 * NY) / 3; ++j) {
        for (size_t i = NX / 3; i < (2 * NX) / 3; ++i) {
            field->values[cell_index(i, j)] = 1.0;
        }
    }
}

static double compute_stable_dt(void) {
    const double dx2 = DX * DX;
    const double dy2 = DY * DY;
    const double denom = 2.0 * DIFFUSIVITY * ((1.0 / dx2) + (1.0 / dy2));

    return SAFETY / denom;
}

static void diffuse_explicit(const scalar_field_t *current, scalar_field_t *next, double dt) {
    const double inv_dx2 = 1.0 / (DX * DX);
    const double inv_dy2 = 1.0 / (DY * DY);

    *next = *current;

    for (size_t j = 1; j < NY - 1; ++j) {
        for (size_t i = 1; i < NX - 1; ++i) {
            const size_t c = cell_index(i, j);
            const double center = current->values[c];
            const double laplacian_x =
                (current->values[cell_index(i + 1, j)] - 2.0 * center + current->values[cell_index(i - 1, j)]) * inv_dx2;
            const double laplacian_y =
                (current->values[cell_index(i, j + 1)] - 2.0 * center + current->values[cell_index(i, j - 1)]) * inv_dy2;

            next->values[c] = center + dt * DIFFUSIVITY * (laplacian_x + laplacian_y);
        }
    }
}

static double l2_change_norm(const scalar_field_t *a, const scalar_field_t *b) {
    double delta[NX * NY];

    for (size_t k = 0; k < NX * NY; ++k) {
        delta[k] = b->values[k] - a->values[k];
    }

    gsl_vector_view delta_view = gsl_vector_view_array(delta, NX * NY);
    return gsl_blas_dnrm2(&delta_view.vector);
}

static double total_scalar(const scalar_field_t *field) {
    double sum = 0.0;

    for (size_t k = 0; k < NX * NY; ++k) {
        sum += field->values[k];
    }

    return sum * DX * DY;
}

static void print_centerline(const scalar_field_t *field) {
    const size_t mid = NY / 2;

    puts("\nCenterline sample (x, phi):");
    for (size_t i = 0; i < NX; i += 5) {
        const double x = i * DX;
        const double phi = field->values[cell_index(i, mid)];
        printf("  %.3f  %.6f\n", x, phi);
    }
}

int main(void) {
    scalar_field_t current;
    scalar_field_t next;
    const double dt = compute_stable_dt();

    field_fill(&current, 0.0);
    apply_boundary_conditions(&current);

    printf("CFD solver starter: 2D diffusion on a %dx%d grid\n", NX, NY);
    printf("dx = %.6f, dy = %.6f, nu = %.6f, dt = %.6f\n", DX, DY, DIFFUSIVITY, dt);

    for (int step = 1; step <= STEPS; ++step) {
        diffuse_explicit(&current, &next, dt);
        apply_boundary_conditions(&next);

        if (step == 1 || step % REPORT_EVERY == 0 || step == STEPS) {
            const double residual = l2_change_norm(&current, &next);
            printf("step %4d  residual(L2) = %.10e  total_scalar = %.10f\n",
                   step,
                   residual,
                   total_scalar(&next));
        }

        current = next;
    }

    print_centerline(&current);

    return EXIT_SUCCESS;
}
