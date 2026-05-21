#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_vector.h>

#include "raylib.h"

#define GRID_FIELD_COUNT 5

typedef struct {
    size_t nx;
    size_t ny;
    double dt;
    double rho;
    double nu;
    double lid_velocity;
    double pressure_tolerance;
    size_t max_pressure_iterations;
    int steps_per_frame;
    int cell_size;
} simulation_config_t;

typedef struct {
    size_t nx_total;
    size_t ny_total;
    size_t cell_count;
    double dx;
    double dy;
} grid_shape_t;

typedef struct {
    double *storage;
    double *p;
    double *u;
    double *v;
    double *us;
    double *vs;
} staggered_grid_t;

typedef struct {
    gsl_spmatrix *laplacian;
    gsl_spmatrix *laplacian_csr;
    gsl_splinalg_itersolve *solver;
    gsl_vector *rhs;
    gsl_vector *pressure;
    size_t last_iterations;
    double last_residual;
} pressure_workspace_t;

typedef struct {
    simulation_config_t config;
    grid_shape_t shape;
    staggered_grid_t grid;
    pressure_workspace_t pressure;
    size_t step_count;
} simulation_t;

static simulation_config_t default_config(void) {
    return (simulation_config_t) {
        .nx = 64,
        .ny = 64,
        .dt = 0.001,
        .rho = 1.0,
        .nu = 0.01,
        .lid_velocity = 1.0,
        .pressure_tolerance = 1e-5,
        .max_pressure_iterations = 250,
        .steps_per_frame = 4,
        .cell_size = 12
    };
}

static void print_usage(const char *program_name) {
    printf("Usage: %s [options]\n", program_name);
    printf("Options:\n");
    printf("  --nx N                 Interior cells in x (default: 64)\n");
    printf("  --ny N                 Interior cells in y (default: 64)\n");
    printf("  --dt VALUE             Time step (default: 0.001)\n");
    printf("  --rho VALUE            Fluid density (default: 1.0)\n");
    printf("  --nu VALUE             Kinematic viscosity (default: 0.01)\n");
    printf("  --lid VALUE            Lid velocity (default: 1.0)\n");
    printf("  --tol VALUE            Pressure solve tolerance (default: 1e-5)\n");
    printf("  --pressure-iters N     Maximum pressure iterations per step (default: 250)\n");
    printf("  --steps N              Solver steps per rendered frame (default: 4)\n");
    printf("  --cell N               Rendered pixel size per cell (default: 12)\n");
    printf("  --help                 Show this message\n");
}

static bool parse_size_option(const char *name, const char *value, size_t min_value, size_t max_value, size_t *out) {
    char *end = NULL;
    errno = 0;
    unsigned long parsed = strtoul(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0' || parsed < min_value || parsed > max_value) {
        fprintf(stderr, "Invalid %s: expected integer in [%zu, %zu], got '%s'.\n", name, min_value, max_value, value);
        return false;
    }

    *out = (size_t) parsed;
    return true;
}

static bool parse_int_option(const char *name, const char *value, int min_value, int max_value, int *out) {
    char *end = NULL;
    errno = 0;
    long parsed = strtol(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0' || parsed < min_value || parsed > max_value) {
        fprintf(stderr, "Invalid %s: expected integer in [%d, %d], got '%s'.\n", name, min_value, max_value, value);
        return false;
    }

    *out = (int) parsed;
    return true;
}

static bool parse_double_option(const char *name, const char *value, double min_value, double max_value, double *out) {
    char *end = NULL;
    errno = 0;
    double parsed = strtod(value, &end);
    if (errno != 0 || end == value || *end != '\0' || !isfinite(parsed) || parsed < min_value || parsed > max_value) {
        fprintf(stderr, "Invalid %s: expected value in [%g, %g], got '%s'.\n", name, min_value, max_value, value);
        return false;
    }

    *out = parsed;
    return true;
}

static bool parse_arguments(int argc, char **argv, simulation_config_t *config) {
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (strcmp(arg, "--help") == 0) {
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
        }

        if (i + 1 >= argc) {
            fprintf(stderr, "Missing value after %s.\n", arg);
            return false;
        }

        const char *value = argv[++i];
        if (strcmp(arg, "--nx") == 0) {
            if (!parse_size_option(arg, value, 4, 512, &config->nx)) {
                return false;
            }
        } else if (strcmp(arg, "--ny") == 0) {
            if (!parse_size_option(arg, value, 4, 512, &config->ny)) {
                return false;
            }
        } else if (strcmp(arg, "--dt") == 0) {
            if (!parse_double_option(arg, value, 1e-8, 0.1, &config->dt)) {
                return false;
            }
        } else if (strcmp(arg, "--rho") == 0) {
            if (!parse_double_option(arg, value, 1e-8, 1000.0, &config->rho)) {
                return false;
            }
        } else if (strcmp(arg, "--nu") == 0) {
            if (!parse_double_option(arg, value, 0.0, 100.0, &config->nu)) {
                return false;
            }
        } else if (strcmp(arg, "--lid") == 0) {
            if (!parse_double_option(arg, value, -100.0, 100.0, &config->lid_velocity)) {
                return false;
            }
        } else if (strcmp(arg, "--tol") == 0) {
            if (!parse_double_option(arg, value, 1e-14, 1.0, &config->pressure_tolerance)) {
                return false;
            }
        } else if (strcmp(arg, "--pressure-iters") == 0) {
            if (!parse_size_option(arg, value, 1, 100000, &config->max_pressure_iterations)) {
                return false;
            }
        } else if (strcmp(arg, "--steps") == 0) {
            if (!parse_int_option(arg, value, 1, 1000, &config->steps_per_frame)) {
                return false;
            }
        } else if (strcmp(arg, "--cell") == 0) {
            if (!parse_int_option(arg, value, 1, 80, &config->cell_size)) {
                return false;
            }
        } else {
            fprintf(stderr, "Unknown option: %s\n", arg);
            print_usage(argv[0]);
            return false;
        }
    }

    return true;
}

static bool checked_mul_size(size_t a, size_t b, size_t *out) {
    if (a != 0 && b > SIZE_MAX / a) {
        return false;
    }

    *out = a * b;
    return true;
}

static bool build_shape(const simulation_config_t *config, grid_shape_t *shape) {
    size_t nx_total = config->nx + 2;
    size_t ny_total = config->ny + 2;
    size_t cell_count = 0;

    if (!checked_mul_size(nx_total, ny_total, &cell_count)) {
        fprintf(stderr, "Grid is too large to allocate safely.\n");
        return false;
    }

    *shape = (grid_shape_t) {
        .nx_total = nx_total,
        .ny_total = ny_total,
        .cell_count = cell_count,
        .dx = 1.0 / (double) config->nx,
        .dy = 1.0 / (double) config->ny
    };

    return true;
}

static size_t ix(const simulation_t *sim, size_t i, size_t j) {
    return j * sim->shape.nx_total + i;
}

static size_t pressure_ix(const simulation_t *sim, size_t i, size_t j) {
    return (j - 1) * sim->config.nx + (i - 1);
}

static bool allocate_grid(const grid_shape_t *shape, staggered_grid_t *grid) {
    size_t field_values = 0;
    if (!checked_mul_size(shape->cell_count, GRID_FIELD_COUNT, &field_values)) {
        fprintf(stderr, "Grid field storage is too large to allocate safely.\n");
        return false;
    }

    grid->storage = calloc(field_values, sizeof(*grid->storage));
    if (grid->storage == NULL) {
        fprintf(stderr, "Failed to allocate %.2f MB for grid storage.\n",
                (double) (field_values * sizeof(*grid->storage)) / (1024.0 * 1024.0));
        return false;
    }

    grid->p = grid->storage;
    grid->u = grid->p + shape->cell_count;
    grid->v = grid->u + shape->cell_count;
    grid->us = grid->v + shape->cell_count;
    grid->vs = grid->us + shape->cell_count;
    return true;
}

static void free_grid(staggered_grid_t *grid) {
    free(grid->storage);
    *grid = (staggered_grid_t) {0};
}

static bool set_matrix_entry(gsl_spmatrix *matrix, size_t row, size_t col, double value) {
    int status = gsl_spmatrix_set(matrix, row, col, value);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to set sparse matrix entry (%zu, %zu): GSL status %d.\n", row, col, status);
        return false;
    }

    return true;
}

static bool build_laplacian(simulation_t *sim) {
    const simulation_config_t *config = &sim->config;
    const double dxi2 = 1.0 / (sim->shape.dx * sim->shape.dx);
    const double dyi2 = 1.0 / (sim->shape.dy * sim->shape.dy);
    size_t pressure_cells = 0;
    size_t max_nonzeros = 0;

    if (!checked_mul_size(config->nx, config->ny, &pressure_cells) ||
        !checked_mul_size(pressure_cells, 5, &max_nonzeros)) {
        fprintf(stderr, "Pressure matrix is too large to allocate safely.\n");
        return false;
    }

    sim->pressure.laplacian = gsl_spmatrix_alloc_nzmax(pressure_cells, pressure_cells, max_nonzeros, GSL_SPMATRIX_TRIPLET);
    if (sim->pressure.laplacian == NULL) {
        fprintf(stderr, "Failed to allocate pressure Laplacian.\n");
        return false;
    }

    for (size_t j = 1; j <= config->ny; ++j) {
        for (size_t i = 1; i <= config->nx; ++i) {
            const size_t row = pressure_ix(sim, i, j);
            double diag = 0.0;

            if (i == 1 && j == 1) {
                if (!set_matrix_entry(sim->pressure.laplacian, row, row, 1.0)) {
                    return false;
                }
                continue;
            }

            if (i > 1) {
                if (!set_matrix_entry(sim->pressure.laplacian, row, pressure_ix(sim, i - 1, j), -dxi2)) {
                    return false;
                }
                diag += dxi2;
            }
            if (i < config->nx) {
                if (!set_matrix_entry(sim->pressure.laplacian, row, pressure_ix(sim, i + 1, j), -dxi2)) {
                    return false;
                }
                diag += dxi2;
            }
            if (j > 1) {
                if (!set_matrix_entry(sim->pressure.laplacian, row, pressure_ix(sim, i, j - 1), -dyi2)) {
                    return false;
                }
                diag += dyi2;
            }
            if (j < config->ny) {
                if (!set_matrix_entry(sim->pressure.laplacian, row, pressure_ix(sim, i, j + 1), -dyi2)) {
                    return false;
                }
                diag += dyi2;
            }
            if (!set_matrix_entry(sim->pressure.laplacian, row, row, diag)) {
                return false;
            }
        }
    }

    sim->pressure.laplacian_csr = gsl_spmatrix_crs(sim->pressure.laplacian);
    if (sim->pressure.laplacian_csr == NULL) {
        fprintf(stderr, "Failed to convert pressure Laplacian to CSR format.\n");
        return false;
    }

    sim->pressure.solver = gsl_splinalg_itersolve_alloc(gsl_splinalg_itersolve_gmres, pressure_cells, 0);
    sim->pressure.rhs = gsl_vector_alloc(pressure_cells);
    sim->pressure.pressure = gsl_vector_alloc(pressure_cells);
    if (sim->pressure.solver == NULL || sim->pressure.rhs == NULL || sim->pressure.pressure == NULL) {
        fprintf(stderr, "Failed to allocate reusable pressure solver workspace.\n");
        return false;
    }

    return true;
}

static void free_pressure_workspace(pressure_workspace_t *workspace) {
    if (workspace->pressure != NULL) {
        gsl_vector_free(workspace->pressure);
    }
    if (workspace->rhs != NULL) {
        gsl_vector_free(workspace->rhs);
    }
    if (workspace->solver != NULL) {
        gsl_splinalg_itersolve_free(workspace->solver);
    }
    if (workspace->laplacian_csr != NULL) {
        gsl_spmatrix_free(workspace->laplacian_csr);
    }
    if (workspace->laplacian != NULL) {
        gsl_spmatrix_free(workspace->laplacian);
    }
    *workspace = (pressure_workspace_t) {0};
}

static bool simulation_init(simulation_t *sim, simulation_config_t config) {
    *sim = (simulation_t) {0};
    sim->config = config;

    if (!build_shape(&sim->config, &sim->shape)) {
        return false;
    }
    if (!allocate_grid(&sim->shape, &sim->grid)) {
        return false;
    }
    if (!build_laplacian(sim)) {
        return false;
    }

    return true;
}

static void simulation_destroy(simulation_t *sim) {
    free_pressure_workspace(&sim->pressure);
    free_grid(&sim->grid);
    *sim = (simulation_t) {0};
}

static void apply_boundary_conditions(simulation_t *sim) {
    staggered_grid_t *grid = &sim->grid;
    const size_t nx = sim->config.nx;
    const size_t ny = sim->config.ny;

    for (size_t j = 0; j < sim->shape.ny_total; ++j) {
        grid->u[ix(sim, 1, j)] = 0.0;
        grid->u[ix(sim, nx + 1, j)] = 0.0;
        grid->v[ix(sim, 0, j)] = -grid->v[ix(sim, 1, j)];
        grid->v[ix(sim, nx + 1, j)] = -grid->v[ix(sim, nx, j)];

        grid->us[ix(sim, 1, j)] = 0.0;
        grid->us[ix(sim, nx + 1, j)] = 0.0;
    }

    for (size_t i = 0; i < sim->shape.nx_total; ++i) {
        grid->v[ix(sim, i, 1)] = 0.0;
        grid->v[ix(sim, i, ny + 1)] = 0.0;
        grid->u[ix(sim, i, 0)] = -grid->u[ix(sim, i, 1)];
        grid->u[ix(sim, i, ny + 1)] = 2.0 * sim->config.lid_velocity - grid->u[ix(sim, i, ny)];

        grid->vs[ix(sim, i, 1)] = 0.0;
        grid->vs[ix(sim, i, ny + 1)] = 0.0;
    }
}

static void compute_tentative_velocity(simulation_t *sim) {
    staggered_grid_t *grid = &sim->grid;
    const simulation_config_t *config = &sim->config;
    const size_t nx = config->nx;
    const size_t ny = config->ny;
    const double dt = config->dt;
    const double nu = config->nu;
    const double dxi = 1.0 / sim->shape.dx;
    const double dyi = 1.0 / sim->shape.dy;
    const double dxi2 = dxi * dxi;
    const double dyi2 = dyi * dyi;

    for (size_t j = 1; j <= ny; ++j) {
        for (size_t i = 2; i <= nx; ++i) {
            const double u_center = grid->u[ix(sim, i, j)];
            const double v_here = 0.25 * (
                grid->v[ix(sim, i - 1, j)] +
                grid->v[ix(sim, i - 1, j + 1)] +
                grid->v[ix(sim, i, j)] +
                grid->v[ix(sim, i, j + 1)]
            );

            const double diff_x = nu * (grid->u[ix(sim, i - 1, j)] - 2.0 * u_center + grid->u[ix(sim, i + 1, j)]) * dxi2;
            const double diff_y = nu * (grid->u[ix(sim, i, j - 1)] - 2.0 * u_center + grid->u[ix(sim, i, j + 1)]) * dyi2;
            const double adv_x = u_center * (grid->u[ix(sim, i + 1, j)] - grid->u[ix(sim, i - 1, j)]) * 0.5 * dxi;
            const double adv_y = v_here * (grid->u[ix(sim, i, j + 1)] - grid->u[ix(sim, i, j - 1)]) * 0.5 * dyi;

            grid->us[ix(sim, i, j)] = u_center + dt * (diff_x + diff_y - adv_x - adv_y);
        }
    }

    for (size_t j = 2; j <= ny; ++j) {
        for (size_t i = 1; i <= nx; ++i) {
            const double v_center = grid->v[ix(sim, i, j)];
            const double u_here = 0.25 * (
                grid->u[ix(sim, i, j - 1)] +
                grid->u[ix(sim, i, j)] +
                grid->u[ix(sim, i + 1, j - 1)] +
                grid->u[ix(sim, i + 1, j)]
            );

            const double diff_x = nu * (grid->v[ix(sim, i - 1, j)] - 2.0 * v_center + grid->v[ix(sim, i + 1, j)]) * dxi2;
            const double diff_y = nu * (grid->v[ix(sim, i, j - 1)] - 2.0 * v_center + grid->v[ix(sim, i, j + 1)]) * dyi2;
            const double adv_x = u_here * (grid->v[ix(sim, i + 1, j)] - grid->v[ix(sim, i - 1, j)]) * 0.5 * dxi;
            const double adv_y = v_center * (grid->v[ix(sim, i, j + 1)] - grid->v[ix(sim, i, j - 1)]) * 0.5 * dyi;

            grid->vs[ix(sim, i, j)] = v_center + dt * (diff_x + diff_y - adv_x - adv_y);
        }
    }
}

static int solve_pressure(simulation_t *sim) {
    staggered_grid_t *grid = &sim->grid;
    pressure_workspace_t *pressure = &sim->pressure;
    const simulation_config_t *config = &sim->config;
    const double rhs_scale = -(config->rho / config->dt);

    for (size_t j = 1; j <= config->ny; ++j) {
        for (size_t i = 1; i <= config->nx; ++i) {
            const size_t row = pressure_ix(sim, i, j);
            const double div =
                (grid->us[ix(sim, i + 1, j)] - grid->us[ix(sim, i, j)]) / sim->shape.dx +
                (grid->vs[ix(sim, i, j + 1)] - grid->vs[ix(sim, i, j)]) / sim->shape.dy;

            gsl_vector_set(pressure->rhs, row, rhs_scale * div);
            gsl_vector_set(pressure->pressure, row, grid->p[ix(sim, i, j)]);
        }
    }

    gsl_vector_set(pressure->rhs, pressure_ix(sim, 1, 1), 0.0);

    int status = GSL_CONTINUE;
    size_t iteration = 0;
    for (; iteration < config->max_pressure_iterations; ++iteration) {
        status = gsl_splinalg_itersolve_iterate(
            pressure->laplacian_csr,
            pressure->rhs,
            config->pressure_tolerance,
            pressure->pressure,
            pressure->solver
        );

        if (status == GSL_SUCCESS) {
            break;
        }
        if (status != GSL_CONTINUE) {
            fprintf(stderr, "Pressure solver failed after %zu iterations: GSL status %d.\n", iteration + 1, status);
            return status;
        }
    }

    pressure->last_iterations = (iteration < config->max_pressure_iterations) ? iteration + 1 : config->max_pressure_iterations;
    pressure->last_residual = gsl_splinalg_itersolve_normr(pressure->solver);

    if (status == GSL_CONTINUE) {
        fprintf(stderr, "Pressure solver reached %zu iterations without meeting tolerance %.3e. Residual: %.3e\n",
                config->max_pressure_iterations,
                config->pressure_tolerance,
                pressure->last_residual);
    }

    for (size_t j = 1; j <= config->ny; ++j) {
        for (size_t i = 1; i <= config->nx; ++i) {
            grid->p[ix(sim, i, j)] = gsl_vector_get(pressure->pressure, pressure_ix(sim, i, j));
        }
    }

    return GSL_SUCCESS;
}

static void apply_corrector(simulation_t *sim) {
    staggered_grid_t *grid = &sim->grid;
    const simulation_config_t *config = &sim->config;
    const double pressure_scale = config->dt / config->rho;
    const double dxi = 1.0 / sim->shape.dx;
    const double dyi = 1.0 / sim->shape.dy;

    for (size_t j = 1; j <= config->ny; ++j) {
        for (size_t i = 2; i <= config->nx; ++i) {
            grid->u[ix(sim, i, j)] =
                grid->us[ix(sim, i, j)] -
                pressure_scale * (grid->p[ix(sim, i, j)] - grid->p[ix(sim, i - 1, j)]) * dxi;
        }
    }

    for (size_t j = 2; j <= config->ny; ++j) {
        for (size_t i = 1; i <= config->nx; ++i) {
            grid->v[ix(sim, i, j)] =
                grid->vs[ix(sim, i, j)] -
                pressure_scale * (grid->p[ix(sim, i, j)] - grid->p[ix(sim, i, j - 1)]) * dyi;
        }
    }
}

static int simulation_step(simulation_t *sim) {
    apply_boundary_conditions(sim);
    compute_tentative_velocity(sim);
    apply_boundary_conditions(sim);

    int pressure_status = solve_pressure(sim);
    if (pressure_status != GSL_SUCCESS) {
        return pressure_status;
    }

    apply_corrector(sim);
    apply_boundary_conditions(sim);
    ++sim->step_count;
    return GSL_SUCCESS;
}

static Color color_for_speed(double speed, double scale) {
    double t = speed / scale;
    if (t < 0.0) {
        t = 0.0;
    } else if (t > 1.0) {
        t = 1.0;
    }

    const unsigned char r = (unsigned char) (255.0 * fmax(0.0, 1.5 * t - 0.25));
    const unsigned char g = (unsigned char) (255.0 * sin(t * 3.14159265358979323846));
    const unsigned char b = (unsigned char) (255.0 * (1.0 - 0.85 * t));
    return (Color) {r, g, b, 255};
}

static void draw_velocity_vectors(const simulation_t *sim) {
    const staggered_grid_t *grid = &sim->grid;
    const int cell_size = sim->config.cell_size;
    const size_t stride = sim->config.nx >= 80 ? 8 : 4;
    const double arrow_scale = 0.35 * (double) cell_size / fmax(fabs(sim->config.lid_velocity), 0.1);

    for (size_t j = stride; j <= sim->config.ny; j += stride) {
        for (size_t i = stride; i <= sim->config.nx; i += stride) {
            const double u_c = 0.5 * (grid->u[ix(sim, i, j)] + grid->u[ix(sim, i + 1, j)]);
            const double v_c = 0.5 * (grid->v[ix(sim, i, j)] + grid->v[ix(sim, i, j + 1)]);
            const float cx = (float) ((i - 0.5) * (double) cell_size);
            const float cy = (float) ((sim->config.ny - j + 0.5) * (double) cell_size);
            const Vector2 start = {cx, cy};
            const Vector2 end = {
                cx + (float) (u_c * arrow_scale),
                cy - (float) (v_c * arrow_scale)
            };

            DrawLineEx(start, end, 1.0f, (Color) {235, 245, 255, 160});
        }
    }
}

static void render_simulation(const simulation_t *sim) {
    const staggered_grid_t *grid = &sim->grid;
    const int cell_size = sim->config.cell_size;
    const double speed_scale = fmax(fabs(sim->config.lid_velocity), 0.2);

    BeginDrawing();
    ClearBackground((Color) {8, 10, 14, 255});

    for (size_t j = 1; j <= sim->config.ny; ++j) {
        for (size_t i = 1; i <= sim->config.nx; ++i) {
            const double u_c = 0.5 * (grid->u[ix(sim, i, j)] + grid->u[ix(sim, i + 1, j)]);
            const double v_c = 0.5 * (grid->v[ix(sim, i, j)] + grid->v[ix(sim, i, j + 1)]);
            const double speed = sqrt(u_c * u_c + v_c * v_c);
            const int x = (int) ((i - 1) * (size_t) cell_size);
            const int y = (int) ((sim->config.ny - j) * (size_t) cell_size);

            DrawRectangle(x, y, cell_size, cell_size, color_for_speed(speed, speed_scale));
        }
    }

    draw_velocity_vectors(sim);

    DrawRectangle(8, 8, 280, 68, (Color) {0, 0, 0, 140});
    DrawText(TextFormat("step: %zu", sim->step_count), 16, 14, 16, RAYWHITE);
    DrawText(TextFormat("residual: %.2e", sim->pressure.last_residual), 16, 34, 16, RAYWHITE);
    DrawText(TextFormat("pressure iterations: %zu", sim->pressure.last_iterations), 16, 54, 16, RAYWHITE);
    DrawFPS(GetScreenWidth() - 90, 12);

    EndDrawing();
}

int main(int argc, char **argv) {
    gsl_set_error_handler_off();

    simulation_config_t config = default_config();
    if (!parse_arguments(argc, argv, &config)) {
        return EXIT_FAILURE;
    }

    simulation_t sim = {0};
    if (!simulation_init(&sim, config)) {
        simulation_destroy(&sim);
        return EXIT_FAILURE;
    }

    const int window_width = (int) (config.nx * (size_t) config.cell_size);
    const int window_height = (int) (config.ny * (size_t) config.cell_size);
    InitWindow(window_width, window_height, "Navier-Stokes Lid-Driven Cavity");
    SetTargetFPS(60);

    int status = EXIT_SUCCESS;
    while (!WindowShouldClose()) {
        for (int step = 0; step < config.steps_per_frame; ++step) {
            if (simulation_step(&sim) != GSL_SUCCESS) {
                status = EXIT_FAILURE;
                break;
            }
        }

        render_simulation(&sim);
        if (status != EXIT_SUCCESS) {
            break;
        }
    }

    CloseWindow();
    simulation_destroy(&sim);
    return status;
}
