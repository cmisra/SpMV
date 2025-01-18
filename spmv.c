#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// Define a structure for COO format
typedef struct {
    int* row;     // Row indices
    int* col;     // Column indices
    double* data; // Non-zero values
    int nnz;      // Number of non-zero elements
    int size;     // Size of the matrix (n x n)
} COOMatrix;

// Define a structure for DIA format
typedef struct {
    int* offsets;  // Diagonal offsets
    double** data; // Diagonal values
    int ndiag;     // Number of diagonals
    int size;      // Size of the matrix (n x n)
} DIAMatrix;

// Function to initialize a COO matrix
void initializeCOOMatrix(COOMatrix* matrix, int size) {
    matrix->size = size;
    matrix->nnz = size; // For a diagonal matrix, nnz equals the size
    matrix->row = (int*)malloc(size * sizeof(int));
    matrix->col = (int*)malloc(size * sizeof(int));
    matrix->data = (double*)malloc(size * sizeof(double));

    if (matrix->row == NULL || matrix->col == NULL || matrix->data == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
}

// Function to initialize a DIA matrix
void initializeDIAMatrix(DIAMatrix* matrix, int size, int ndiag) {
    matrix->size = size;
    matrix->ndiag = ndiag;
    matrix->offsets = (int*)malloc(ndiag * sizeof(int));
    matrix->data = (double**)malloc(ndiag * sizeof(double*));
    for (int i = 0; i < ndiag; i++) {
        matrix->data[i] = (double*)calloc(size, sizeof(double));
    }

    if (matrix->offsets == NULL || matrix->data == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
}

// Function to free the memory allocated for the DIA matrix
void freeDIAMatrix(DIAMatrix* matrix) {
    free(matrix->offsets);
    for (int i = 0; i < matrix->ndiag; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
}

// Function to convert COO to DIA format
void cooToDia(const COOMatrix* coo, DIAMatrix* dia) {
    int max_diags = 2 * coo->size - 1; // Maximum number of diagonals
    int* diag_counts = (int*)calloc(max_diags, sizeof(int));
    int base_offset = coo->size - 1;

    // Count diagonals
    for (int i = 0; i < coo->nnz; i++) {
        int diag = coo->col[i] - coo->row[i] + base_offset;
        diag_counts[diag]++;
    }

    // Determine the number of active diagonals
    int ndiag = 0;
    for (int i = 0; i < max_diags; i++) {
        if (diag_counts[i] > 0) {
            ndiag++;
        }
    }

    // Initialize DIA matrix
    initializeDIAMatrix(dia, coo->size, ndiag);

    // Store active diagonals and their offsets
    int diag_index = 0;
    for (int i = 0; i < max_diags; i++) {
        if (diag_counts[i] > 0) {
            dia->offsets[diag_index] = i - base_offset;
            diag_index++;
        }
    }

    // Populate DIA matrix with data
    for (int i = 0; i < coo->nnz; i++) {
        int row = coo->row[i];
        int col = coo->col[i];
        int diag = col - row + base_offset;

        for (int j = 0; j < dia->ndiag; j++) {
            if (dia->offsets[j] == diag - base_offset) {
                dia->data[j][row] = coo->data[i];
                break;
            }
        }
    }

    free(diag_counts);
}

// Function to perform Sparse Matrix-Vector Multiplication (SpMV) using DIA format
void spmvDia(const DIAMatrix* matrix, const double* vector, double* result) {
    // Initialize the result vector to zero
    for (int i = 0; i < matrix->size; i++) {
        result[i] = 0.0;
    }

    // Perform SpMV
    for (int d = 0; d < matrix->ndiag; d++) {
        int offset = matrix->offsets[d];
        for (int i = 0; i < matrix->size; i++) {
            int j = i + offset;
            if (j >= 0 && j < matrix->size) {
                result[i] += matrix->data[d][i] * vector[j];
            }
        }
    }
}

// Function to generate a random vector
void generateRandomVector(double* vector, int size) {
    for (int i = 0; i < size; i++) {
        vector[i] = (rand() % 10000) / 100.0; // Random double value between 0.0 and 99.99
    }
}

// Function to generate a random diagonal matrix in COO format
void generateDiagonalMatrix(COOMatrix* matrix) {
    for (int i = 0; i < matrix->size; i++) {
        matrix->row[i] = i;          // Row index (diagonal elements)
        matrix->col[i] = i;          // Column index
        matrix->data[i] = (rand() % 10000) / 100.0; // Random double value between 0.00 and 99.99
    }
}

// Function to print the COO matrix
void printCOOMatrix(const COOMatrix* matrix) {
    printf("COO Format:\n");
    printf("Row\tCol\tValue\n");
    for (int i = 0; i < matrix->nnz; i++) {
        printf("%d\t%d\t%.2f\n", matrix->row[i], matrix->col[i], matrix->data[i]);
    }
}

// Function to free the memory allocated for the COO matrix
void freeCOOMatrix(COOMatrix* matrix) {
    free(matrix->row);
    free(matrix->col);
    free(matrix->data);
}

// Function to perform Sparse Matrix-Vector Multiplication (SpMV)
void spmv(const COOMatrix* matrix, const double* vector, double* result) { 

    // Initialize the result vector to zero
    for (int i = 0; i < matrix->size; i++) {
        result[i] = 0.0;
    }

    // Perform SpMV: result[row[i]] += data[i] * vector[col[i]]
    for (int i = 0; i < matrix->nnz; i++) {
        result[matrix->row[i]] += matrix->data[i] * vector[matrix->col[i]];
    }
}

// Function to perform Sparse Matrix-Vector Multiplication (SpMV) using COO format and OpenMP
void spmvCooOpenMP(const COOMatrix* matrix, const double* vector, double* result) {
    // Set the number of threads to 4
    omp_set_num_threads(4);

    // Initialize the result vector to zero
    int i = 0;
    #pragma omp parallel for
    for (i = 0; i < matrix->size; i++) {
        result[i] = 0.0;
    }

    // Perform SpMV: result[row[i]] += data[i] * vector[col[i]]
    #pragma omp parallel for schedule(static)
    for (i = 0; i < matrix->nnz; i++) {
        #pragma omp atomic
        result[matrix->row[i]] += matrix->data[i] * vector[matrix->col[i]];
    }
}

int main() {
    int size;            // Size of the square matrix
    COOMatrix cooMatrix; // COO matrix structure
    DIAMatrix diaMatrix;
    double* vector, * result;
    clock_t start, end; // Variables for timing

    // Seed the random number generator
    srand(time(NULL));

    // Input matrix size
    printf("Enter the size of the diagonal matrix (n x n): ");
    scanf("%d", &size);

    // Initialize and generate the COO matrix
    initializeCOOMatrix(&cooMatrix, size);
    generateDiagonalMatrix(&cooMatrix);

    // Convert COO to DIA
    //cooToDia(&cooMatrix, &diaMatrix);

    // Allocate memory for the dense vector and result
    vector = (double*)malloc(size * sizeof(double));
    result = (double*)malloc(size * sizeof(double));

    if (vector == NULL || result == NULL) {
        printf("Memory allocation failed.\n");
        freeCOOMatrix(&cooMatrix);
        exit(1);
    }

    // Generate a random vector
    generateRandomVector(vector, size);

    // Measure the time taken for SpMV
    /*start = clock();
    spmv(&cooMatrix, vector, result);
    end = clock();*/

    // Calculate and display the time taken
    /*double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for SpMV: %.6f seconds\n", time_taken);*/

    // Measure the time taken for SpMV
    start = clock();
    spmvCooOpenMP(&cooMatrix, vector, result);
    end = clock();

    // Calculate and display the time taken
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for SpMV in (COO and OpenMP): %.6f seconds\n", time_taken);

    // Free the allocated memory
    freeCOOMatrix(&cooMatrix);

    // Measure the time taken for SpMV using DIA format
    /*start = clock();
    spmvDia(&diaMatrix, vector, result);
    end = clock();*/

    // Calculate and display the time taken
    /*double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for SpMV (DIA): %.6f seconds\n", time_taken);*/

    // Free the allocated memory
    //freeDIAMatrix(&diaMatrix);
    free(vector);
    free(result);

    return 0;
}
