#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



double** matrix_unconst( double theta[], int rows, int cols) {
  double** matrix = (double**)malloc(rows * sizeof(double*));
  for (int i = 0; i < rows; i++) {
    matrix[i] = (double*)malloc(cols * sizeof(double));
  }

  int k=0;
  for( int i=0; i<rows; i++){
    for(int j=0; j<cols;j++){
      matrix[i][j]=theta[k];
      k=k+1;
    }
  }

  return matrix;
}




int* find_non_zero_indices(double *arr, int d, int *num_non_zero) {
    // Allocate memory to store indices (maximum size is d)
    int *indices = (int *)malloc(d * sizeof(int));
    *num_non_zero = 0;  // Initialize the count of non-zero elements

    // Loop through the array to find non-zero elements
    for (int i = 0; i < d; i++) {

        if (arr[i] != 0) {
            indices[*num_non_zero] = i;  // Store the index of the non-zero element
            (*num_non_zero)++;           // Increment the count
        }
    }
    return indices;
}

double pnorm(double x) {
    return 0.5 * (1 + erf(x / sqrt(2.0)));
}


double bi_stdf_HR(double * x , double Gamma){
    int num_non_zero = 0;
    int *indices = find_non_zero_indices(x, 2, &num_non_zero);
  if(num_non_zero== 1) {
    return(x[indices[0]  ]);
  }
  if(num_non_zero== 0) {
    return(0.0);
  }
    double term1 = x[0]*pnorm( (log(x[0] / x[1] )/ sqrt(Gamma)) +  sqrt(Gamma)/2.0 ) ;
    double term2 = x[1]*pnorm( (log(x[1] / x[0] )/ sqrt(Gamma)) +  sqrt(Gamma)/2.0 );
    return term1 + term2;
}




double stdf_log(int d, int r, double  **A, double alpha[r], double x[d]){
  double xA[d][r];
  for (int i=0;i <d; i++){
    for (int j=0; j<r; j++){
      xA[i][j]=pow(x[i]* A[i][j], 1/alpha[j]);
    }
  }
  double scol;
  double s=0;
  for( int j=0; j<r; j++){
    scol=0;
    for (int i=0; i<d; i++){
      scol=scol+xA[i][j];
    }
    s=s+pow(scol, alpha[j]);
  }
  return s;
}



double norm_p_element(int d, double vec[d], double p){
double x=0;
for ( int i=0; i<d;i++){
    x=x+pow(fabs(vec[i]),p );
}
x=pow(x, 1/p);
return x;
}


double norm_p_element_2(int d, double vec[d], double p){
double x=0;
for ( int i=0; i<d;i++){
    x=x+pow(fabs(vec[i]),p );
}
return x;
}

double norm_p_matrix(int d, int r, double matrix[d][r], double p){
double x=0;
for ( int i=0; i<d;i++){
for (int j=0; j<r;j++){
    x=x+pow(fabs(matrix[i][j]),p );
}
}
x=pow(x, 1/p);
return x;
}



void SSR_col(double* p,double * lambda, double* theta, int* d, int* k, int* q, double* alpha, double* w, double* Grid_points, double *R) {
  double **M = matrix_unconst(theta, *d, *k);
  double x[*d];
  double row_vec[*d];
  *R=0;
   double flatM[*d][*k];
    for (int i = 0; i < *d; i++) {
        for (int j = 0; j < *k; j++) {
            flatM[i][j] = M[i][j];
        }
    }
  for (int m = 0; m < *q; m++) {
    for( int i=0; i<*d; i++){
      x[i]=Grid_points[(*d)* (m)+i];
    }
    *R = *R + pow(w[m] - stdf_log(*d,*k,M,alpha,x), 2);
  }
     for( int j=0; j<*k; j++){
    for(int i=0; i<*d; i++){
        row_vec[i]=flatM[i][j];
    }
   *R=*R+ (*lambda)* (norm_p_element(*d, row_vec,*p));
  }
}


void SSR_matrix(double *p,double *lambda, double* theta, int* d, int* k, int* q, double* alpha, double* w, double* Grid_points, double *R) {
  double **M = matrix_unconst(theta, *d, *k);
  double x[*d];
  *R=0;
  double flatM[*d][*k];
    for (int i = 0; i < *d; i++) {
        for (int j = 0; j < *k; j++) {
            flatM[i][j] = M[i][j];
        }
    }
  for (int m = 0; m < *q; m++) {
    for( int i=0; i<*d; i++){
      x[i]=Grid_points[(*d)* (m)+i];
    }
    *R = *R + pow(w[m] - stdf_log(*d,*k,M,alpha,x), 2);
  }
*R=*R+(*lambda)* norm_p_matrix(*d,*k, flatM,*p);
}









void print_vector(int n, double *vector) {
    for (int i = 0; i < n; i++) {
        printf("%f ", vector[i]);
    }
    printf("\n");  // Newline after printing the vector
}

void print_matrix(int rows, int cols, double **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}






void SSR_row_log(double* p,double * lambda, double* theta, int* d, int* k, int* q, double* alpha, double* w, double* Grid_points, double *R) {
  double **M = matrix_unconst(theta, *d, *k);
  double x[*d];
  *R=0;
  double flatM[*d][*k];
    for (int i = 0; i < *d; i++) {
        for (int j = 0; j < *k; j++) {
            flatM[i][j] = M[i][j];
        }
    }
  for (int m = 0; m < *q; m++) {
    for( int i=0; i<*d; i++){
      x[i]=Grid_points[(*d)* (m)+i];
    }
    *R = *R + pow(w[m] - stdf_log(*d,*k,M,alpha,x), 2);
  }
     for( int i=0; i<*d; i++){
   *R=*R+ (*lambda)* (norm_p_element(*k, flatM[i],*p));
  }
}


double bi_stdf_mix_HR(int d, double x[d], int r, double A[d][r], double ** Gamma) {
    double R = 0.0;
    int num_non_zero = 0;
    int *indices = find_non_zero_indices(x, d, &num_non_zero);

    if (num_non_zero == 0) {
        free(indices);
        return R;
    }

    // We expect bi_stdf_HR to be called only with exactly 2 elements
    if (num_non_zero != 2) {
        free(indices);
        // You can handle this case differently if needed
        return 0.0; // or some safe default
    }

    double xA[num_non_zero][r];
    int row_index;

    for (int i = 0; i < num_non_zero; i++) {
        row_index = indices[i];
        for (int j = 0; j < r; j++) {
            xA[i][j] = x[row_index] * A[row_index][j];
        }
    }

    for (int j = 0; j < r; j++) {
        double vec[2] = {0.0, 0.0};
        for (int i = 0; i < 2; i++) {
            vec[i] = xA[i][j];
        }
        R += bi_stdf_HR(vec, Gamma[indices[0]][indices[1]]);
    }

    free(indices);
    return R;
}


double **construct_matrix(int d, double *vec) {
    // Allocate memory for the matrix (array of row pointers)
    double **matrix = malloc(d * sizeof(double *));
    if (!matrix) {
        fprintf(stderr, "Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row and populate the matrix
    for (int i = 0; i < d; i++) {
        matrix[i] = malloc(d * sizeof(double));
        if (!matrix[i]) {
            fprintf(stderr, "Memory allocation error.\n");
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < d; j++) {
            matrix[i][j] = vec[i * d + j];
        }
    }

    return matrix;
}



void free_flat_matrix(double **mat) {
    if (mat != NULL) {
        free(mat[0]);  // Frees the contiguous data block
        free(mat);     // Frees the array of row pointers
    }
}


double **construct_symmetric_matrix(double *vec, int vec_length) {
    int d = (1 + sqrt(1 + 8 * vec_length)) / 2;
    if (d * (d - 1) / 2 != vec_length) {
        fprintf(stderr, "Vector length is incorrect for a symmetric matrix.\n");
        exit(EXIT_FAILURE);
    }

    double *data = malloc(d * d * sizeof(double));
    if (!data) {
        fprintf(stderr, "Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }

    double **mat = malloc(d * sizeof(double *));
    if (!mat) {
        free(data);
        fprintf(stderr, "Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < d; i++) {
        mat[i] = data + i * d;
        mat[i][i] = 1.0;
    }

    int index = 0;
    for (int i = 0; i < d - 1; i++) {
        for (int j = i + 1; j < d; j++) {
            mat[i][j] = vec[index++];
        }
    }

    for (int i = 0; i < d; i++) {
        for (int j = 0; j < i; j++) {
            mat[i][j] = mat[j][i];
        }
    }

    return mat;
}







void free_matrix(double **mat, int d) {
    for (int i = 0; i < d; i++) {
        free(mat[i]);
    }
    free(mat);
}


void SSR_row_HR(double* p, double* lambda, double* theta, int* d, int* k, int* q,
                double* Gamma, double* w, double* Grid_points, double* R) {
    
    int vec_length = (*d) * ((*d) - 1) / 2;

    // Construct symmetric matrix from vector
    double** Gamma_mat = construct_symmetric_matrix(Gamma, vec_length);

    // Convert theta to matrix form
    double** M = matrix_unconst(theta, *d, *k);

    *R = 0.0;

    // Allocate contiguous memory for flatM
    double (*flatM)[*k] = malloc(*d * sizeof(*flatM));
    if (!flatM) {
        fprintf(stderr, "Memory allocation error for flatM.\n");
        free_flat_matrix(Gamma_mat);
        free_matrix(M, *d);
        exit(EXIT_FAILURE);
    }

    // Copy M into flatM
    for (int i = 0; i < *d; i++) {
        for (int j = 0; j < *k; j++) {
            flatM[i][j] = M[i][j];
        }
    }

    // Evaluate sum of squared residuals over all grid points
    double x[*d];
    for (int m = 0; m < *q; m++) {
        for (int i = 0; i < *d; i++) {
            x[i] = Grid_points[(*d) * m + i];
        }

        double interm = bi_stdf_mix_HR(*d, x, *k, flatM, Gamma_mat);
        *R += pow(w[m] - interm, 2);
    }

    // Regularization term
    for (int i = 0; i < *d; i++) {
        *R += (*lambda) * norm_p_element(*k, flatM[i], *p);
    }

    // Free memory
    free(flatM);
    free_flat_matrix(Gamma_mat);
    free_matrix(M, *d);
}

//int main () {
//  double p = 0.2;
//  double lambda = 0;
//  double theta[16] = {1.0/3.0 , 0 , 1.0/3.0 , 1.0/3.0 , 0 , 0.5 , 0.5 , 0 , 3.0/4.0 , 0 , 0 , 1.0/4.0 , 1.0/2.0 , 1.0/2.0 , 0 , 0 };
//  int d = 4;
 // int k = 4;
//  int q = 1;
//  double Gamma[16] = { 0.3 , 0.3 , 0.3 , 0.3  , 0.3 , 0.3 };
//  double w = 0.5;
//  double Grid_points[4] = {1 , 0.5 , 0 , 0};
 // double R = 0;
//  SSR_row_HR(&p,&lambda, theta, &d, &k, &q, Gamma, &w, Grid_points, &R);
 // printf("result= %f", R);
//}
