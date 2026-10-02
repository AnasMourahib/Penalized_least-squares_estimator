#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================== Matrix helpers ============================== */

double **matrix_unconst(double theta[], int rows, int cols) {
    double **matrix = (double **) malloc((size_t) rows * sizeof(double *));
    if (matrix == NULL) error("Memory allocation failed for matrix.");

    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *) malloc((size_t) cols * sizeof(double));
        if (matrix[i] == NULL) {
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            error("Memory allocation failed for matrix row.");
        }
    }

    int index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) matrix[i][j] = theta[index++];
    }
    return matrix;
}

void free_matrix(double **mat, int rows) {
    if (mat == NULL) return;
    for (int i = 0; i < rows; i++) free(mat[i]);
    free(mat);
}

void free_flat_matrix(double **mat) {
    if (mat != NULL) {
        free(mat[0]);
        free(mat);
    }
}

int *find_non_zero_indices(double *arr, int d, int *num_non_zero) {
    int *indices = (int *) malloc((size_t) d * sizeof(int));
    if (indices == NULL) error("Memory allocation failed for indices.");
    *num_non_zero = 0;
    for (int i = 0; i < d; i++) {
        if (arr[i] != 0.0) indices[(*num_non_zero)++] = i;
    }
    return indices;
}

double **construct_matrix(int d, double *vec) {
    double **matrix = (double **) malloc((size_t) d * sizeof(double *));
    if (matrix == NULL) error("Memory allocation failed for matrix.");

    for (int i = 0; i < d; i++) {
        matrix[i] = (double *) malloc((size_t) d * sizeof(double));
        if (matrix[i] == NULL) {
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            error("Memory allocation failed for matrix row.");
        }
        for (int j = 0; j < d; j++) matrix[i][j] = vec[i * d + j];
    }
    return matrix;
}

double **construct_symmetric_matrix(double *vec, int vec_length) {
    int d = (int) ((1.0 + sqrt(1.0 + 8.0 * vec_length)) / 2.0);
    if (d * (d - 1) / 2 != vec_length) {
        error("Vector length is incorrect for a symmetric matrix.");
    }

    double *data = (double *) calloc((size_t) d * d, sizeof(double));
    if (data == NULL) error("Memory allocation failed for symmetric matrix.");

    double **mat = (double **) malloc((size_t) d * sizeof(double *));
    if (mat == NULL) {
        free(data);
        error("Memory allocation failed for matrix pointers.");
    }

    for (int i = 0; i < d; i++) {
        mat[i] = data + i * d;
        mat[i][i] = 0.0;
    }

    int index = 0;
    for (int i = 0; i < d - 1; i++) {
        for (int j = i + 1; j < d; j++) mat[i][j] = vec[index++];
    }
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < i; j++) mat[i][j] = mat[j][i];
    }
    return mat;
}

/* ==================== Gaussian probability via mvtnorm ===================== */

static void mvtnorm_C_mvtdst(
    int *n, int *nu, double *lower, double *upper,
    int *infin, double *correl, double *delta,
    int *maxpts, double *abseps, double *releps,
    double *error_value, double *value, int *inform, int *rnd
) {
    typedef void (*mvtdst_ptr)(
        int *, int *, double *, double *, int *, double *, double *,
        int *, double *, double *, double *, double *, int *, int *
    );
    static mvtdst_ptr fun = NULL;
    if (fun == NULL) fun = (mvtdst_ptr) R_GetCCallable("mvtnorm", "C_mvtdst");
    fun(n, nu, lower, upper, infin, correl, delta, maxpts, abseps,
        releps, error_value, value, inform, rnd);
}

double mvnorm_cdf_genz(int m, const double *upper, const double *Sigma) {
    if (m == 0) return 1.0;
    if (m == 1) {
        double variance = Sigma[0];
        if (variance <= 0.0 || !R_FINITE(variance)) return NA_REAL;
        return pnorm5(upper[0] / sqrt(variance), 0.0, 1.0, 1, 0);
    }

    double *lower = (double *) calloc((size_t) m, sizeof(double));
    double *upper_std = (double *) calloc((size_t) m, sizeof(double));
    double *delta = (double *) calloc((size_t) m, sizeof(double));
    double *sd = (double *) calloc((size_t) m, sizeof(double));
    int *infin = (int *) calloc((size_t) m, sizeof(int));
    int nc = m * (m - 1) / 2;
    double *correl = (double *) calloc((size_t) nc, sizeof(double));

    if (!lower || !upper_std || !delta || !sd || !infin || !correl) {
        free(lower); free(upper_std); free(delta); free(sd); free(infin); free(correl);
        error("Memory allocation failed in mvnorm_cdf_genz.");
    }

    for (int i = 0; i < m; i++) {
        double variance = Sigma[i * m + i];
        if (variance <= 0.0 || !R_FINITE(variance)) {
            free(lower); free(upper_std); free(delta); free(sd); free(infin); free(correl);
            return NA_REAL;
        }
        sd[i] = sqrt(variance);
        upper_std[i] = upper[i] / sd[i];
        infin[i] = 0;
    }

    int index = 0;
    for (int column = 0; column < m - 1; column++) {
        for (int row = column + 1; row < m; row++) {
            double rho = Sigma[row * m + column] / (sd[row] * sd[column]);
            if (rho > 1.0) rho = 1.0;
            if (rho < -1.0) rho = -1.0;
            correl[index++] = rho;
        }
    }

    int nu = 0, maxpts = 25000 * m, inform = 0, rnd = 1;
    double abseps = 1e-8, releps = 1e-8, estimated_error = 0.0, value = 0.0;

    mvtnorm_C_mvtdst(&m, &nu, lower, upper_std, infin, correl, delta,
                     &maxpts, &abseps, &releps, &estimated_error,
                     &value, &inform, &rnd);

    free(lower); free(upper_std); free(delta); free(sd); free(infin); free(correl);
    if (inform != 0) warning("mvtdst accuracy warning; inform = %d", inform);
    return value;
}

/* ==================== Stable tail dependence functions ===================== */

double normal_cdf(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

double bi_stdf_HR(double *x, double Gamma) {
    int num_non_zero = 0;
    int *indices = find_non_zero_indices(x, 2, &num_non_zero);
    if (num_non_zero == 0) { free(indices); return 0.0; }
    if (num_non_zero == 1) {
        double result = x[indices[0]];
        free(indices);
        return result;
    }
    double sg = sqrt(Gamma);
    double term1 = x[0] * normal_cdf(log(x[0] / x[1]) / sg + sg / 2.0);
    double term2 = x[1] * normal_cdf(log(x[1] / x[0]) / sg + sg / 2.0);
    free(indices);
    return term1 + term2;
}

double stdf_log(int d, int r, double **A, double alpha[r], double x[d]) {
    double xA[d][r];
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < r; j++) xA[i][j] = pow(x[i] * A[i][j], 1.0 / alpha[j]);
    }
    double result = 0.0;
    for (int j = 0; j < r; j++) {
        double s = 0.0;
        for (int i = 0; i < d; i++) s += xA[i][j];
        result += pow(s, alpha[j]);
    }
    return result;
}

double stdf_HR_d(int d, const double *y, double **Gamma) {
    if (d < 1) return 0.0;
    double result = 0.0;

    for (int j = 0; j < d; j++) {
        if (y[j] <= 0.0) continue;
        int m = d - 1;
        if (m == 0) { result += y[j]; continue; }

        double *eta = (double *) calloc((size_t) m, sizeof(double));
        double *Sigma = (double *) calloc((size_t) m * m, sizeof(double));
        int *other = (int *) calloc((size_t) m, sizeof(int));
        if (!eta || !Sigma || !other) {
            free(eta); free(Sigma); free(other);
            error("Memory allocation failed in stdf_HR_d.");
        }

        int pos = 0;
        for (int s = 0; s < d; s++) if (s != j) other[pos++] = s;
        for (int a = 0; a < m; a++) {
            int s = other[a];
            eta[a] = y[s] <= 0.0 ? R_PosInf : log(y[j] / y[s]) + Gamma[j][s] / 2.0;
        }
        for (int a = 0; a < m; a++) {
            int s = other[a];
            for (int b = 0; b < m; b++) {
                int t = other[b];
                Sigma[a * m + b] = (Gamma[j][s] + Gamma[j][t] - Gamma[s][t]) / 2.0;
            }
        }

        double cdf = mvnorm_cdf_genz(m, eta, Sigma);
        free(eta); free(Sigma); free(other);
        if (!R_FINITE(cdf)) return NA_REAL;
        result += y[j] * cdf;
    }
    return result;
}

double stdf_mix_HR_d(int d, const double *x, int r, double A[d][r], double **Gamma) {
    double result = 0.0;
    for (int k = 0; k < r; k++) {
        int size_J = 0;
        for (int j = 0; j < d; j++) if (A[j][k] > 0.0) size_J++;
        if (size_J == 0) continue;

        int *J = (int *) calloc((size_t) size_J, sizeof(int));
        double *y_J = (double *) calloc((size_t) size_J, sizeof(double));
        double **Gamma_J = (double **) calloc((size_t) size_J, sizeof(double *));
        if (!J || !y_J || !Gamma_J) {
            free(J); free(y_J); free(Gamma_J);
            error("Memory allocation failed in stdf_mix_HR_d.");
        }
        for (int s = 0; s < size_J; s++) {
            Gamma_J[s] = (double *) calloc((size_t) size_J, sizeof(double));
            if (!Gamma_J[s]) {
                for (int t = 0; t < s; t++) free(Gamma_J[t]);
                free(Gamma_J); free(y_J); free(J);
                error("Memory allocation failed for Gamma_J.");
            }
        }

        int pos = 0;
        for (int j = 0; j < d; j++) if (A[j][k] > 0.0) J[pos++] = j;
        for (int s = 0; s < size_J; s++) y_J[s] = A[J[s]][k] * x[J[s]];
        for (int s = 0; s < size_J; s++) {
            for (int t = 0; t < size_J; t++) Gamma_J[s][t] = Gamma[J[s]][J[t]];
        }

        double component = stdf_HR_d(size_J, y_J, Gamma_J);
        for (int s = 0; s < size_J; s++) free(Gamma_J[s]);
        free(Gamma_J); free(y_J); free(J);
        if (!R_FINITE(component)) return NA_REAL;
        result += component;
    }
    return result;
}

double bi_stdf_mix_HR(int d, double x[d], int r, double A[d][r], double **Gamma) {
    double result = 0.0;
    int nnz = 0;
    int *indices = find_non_zero_indices(x, d, &nnz);
    if (nnz == 0) { free(indices); return 0.0; }
    if (nnz != 2) { free(indices); return 0.0; }
    for (int k = 0; k < r; k++) {
        double vec[2] = {x[indices[0]] * A[indices[0]][k],
                         x[indices[1]] * A[indices[1]][k]};
        result += bi_stdf_HR(vec, Gamma[indices[0]][indices[1]]);
    }
    free(indices);
    return result;
}

/* Wrapper callable from R through .C(). A_flat must be row-major. */
void test_stdf_mix_HR_d(double *x, int *d, int *r, double *A_flat,
                        double *Gamma_vector, double *result) {
    double **Ap = matrix_unconst(A_flat, *d, *r);
    double (*A)[*r] = (double (*)[*r]) malloc((size_t) *d * sizeof(*A));
    if (!A) { free_matrix(Ap, *d); error("Memory allocation failed for A."); }
    for (int i = 0; i < *d; i++) for (int k = 0; k < *r; k++) A[i][k] = Ap[i][k];

    int glen = (*d) * ((*d) - 1) / 2;
    double **Gamma = construct_symmetric_matrix(Gamma_vector, glen);
    *result = stdf_mix_HR_d(*d, x, *r, A, Gamma);

    free(A);
    free_matrix(Ap, *d);
    free_flat_matrix(Gamma);
}

/* ================================= Norms =================================== */

double norm_p_element(int d, double vec[d], double p) {
    double result = 0.0;
    for (int i = 0; i < d; i++) result += pow(fabs(vec[i]), p);
    return pow(result, 1.0 / p);
}

double norm_p_element_2(int d, double vec[d], double p) {
    double result = 0.0;
    for (int i = 0; i < d; i++) result += pow(fabs(vec[i]), p);
    return result;
}

double norm_p_matrix(int d, int r, double matrix[d][r], double p) {
    double result = 0.0;
    for (int i = 0; i < d; i++) for (int j = 0; j < r; j++) result += pow(fabs(matrix[i][j]), p);
    return pow(result, 1.0 / p);
}

/* =========================== Penalized SSR functions ======================= */

void SSR_col(double *p, double *lambda, double *theta, int *d, int *k, int *q,
             double *alpha, double *w, double *Grid_points, double *R) {
    double **M = matrix_unconst(theta, *d, *k);
    double x[*d], row_vec[*d], flatM[*d][*k];
    *R = 0.0;
    for (int i = 0; i < *d; i++) for (int j = 0; j < *k; j++) flatM[i][j] = M[i][j];
    for (int m = 0; m < *q; m++) {
        for (int i = 0; i < *d; i++) x[i] = Grid_points[(*d) * m + i];
        *R += pow(w[m] - stdf_log(*d, *k, M, alpha, x), 2.0);
    }
    for (int j = 0; j < *k; j++) {
        for (int i = 0; i < *d; i++) row_vec[i] = flatM[i][j];
        *R += (*lambda) * norm_p_element(*d, row_vec, *p);
    }
    free_matrix(M, *d);
}

void SSR_matrix(double *p, double *lambda, double *theta, int *d, int *k, int *q,
                double *alpha, double *w, double *Grid_points, double *R) {
    double **M = matrix_unconst(theta, *d, *k);
    double x[*d], flatM[*d][*k];
    *R = 0.0;
    for (int i = 0; i < *d; i++) for (int j = 0; j < *k; j++) flatM[i][j] = M[i][j];
    for (int m = 0; m < *q; m++) {
        for (int i = 0; i < *d; i++) x[i] = Grid_points[(*d) * m + i];
        *R += pow(w[m] - stdf_log(*d, *k, M, alpha, x), 2.0);
    }
    *R += (*lambda) * norm_p_matrix(*d, *k, flatM, *p);
    free_matrix(M, *d);
}

void SSR_row_log(double *p, double *lambda, double *theta, int *d, int *k, int *q,
                 double *alpha, double *w, double *Grid_points, double *R) {
    double **M = matrix_unconst(theta, *d, *k);
    double x[*d], flatM[*d][*k];
    *R = 0.0;
    for (int i = 0; i < *d; i++) for (int j = 0; j < *k; j++) flatM[i][j] = M[i][j];
    for (int m = 0; m < *q; m++) {
        for (int i = 0; i < *d; i++) x[i] = Grid_points[(*d) * m + i];
        *R += pow(w[m] - stdf_log(*d, *k, M, alpha, x), 2.0);
    }
    for (int i = 0; i < *d; i++) *R += (*lambda) * norm_p_element(*k, flatM[i], *p);
    free_matrix(M, *d);
}

void SSR_row_HR(double *p, double *lambda, double *theta, int *d, int *k, int *q,
                double *Gamma, double *w, double *Grid_points, double *R) {
    int glen = (*d) * ((*d) - 1) / 2;
    double **Gamma_mat = construct_symmetric_matrix(Gamma, glen);
    double **M = matrix_unconst(theta, *d, *k);
    double (*flatM)[*k] = (double (*)[*k]) malloc((size_t) *d * sizeof(*flatM));
    if (!flatM) {
        free_flat_matrix(Gamma_mat); free_matrix(M, *d);
        error("Memory allocation failed for flatM.");
    }
    for (int i = 0; i < *d; i++) for (int j = 0; j < *k; j++) flatM[i][j] = M[i][j];

    *R = 0.0;
    double x[*d];
    for (int m = 0; m < *q; m++) {
        for (int i = 0; i < *d; i++) x[i] = Grid_points[(*d) * m + i];
        double interm = stdf_mix_HR_d(*d, x, *k, flatM, Gamma_mat);
        if (!R_FINITE(interm)) { *R = 1e16; break; }
        *R += pow(w[m] - interm, 2.0);
    }
    for (int i = 0; i < *d; i++) *R += (*lambda) * norm_p_element(*k, flatM[i], *p);

    free(flatM);
    free_flat_matrix(Gamma_mat);
    free_matrix(M, *d);
}

void print_vector(int n, double *vector) {
    for (int i = 0; i < n; i++) printf("%f ", vector[i]);
    printf("\n");
}

void print_matrix(int rows, int cols, double **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) printf("%f ", matrix[i][j]);
        printf("\n");
    }
}
