#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

// Structure contenant les vecteurs du solveur et le stockage CSR de la matrice
typedef struct {
    // Vecteurs pour l'algorithme itératif
    double *R;      // Résidu
    double *D;      // Direction de descente
    double *S;      // Stockage temporaire ou pour d'autres vecteurs
    double *X;      // Solution
    double error;   // Erreur courante
    int size;       // Taille du système (nombre de degrés de liberté)
    int iter;       // Nombre d'itérations effectuées

    // Stockage CSR de la matrice du système
    double *values;   // Tableau des valeurs non nulles
    int *col_index;   // Indices de colonnes associés à chaque valeur
    int *row_ptr;     // Indices de début de chaque ligne (taille = nrows + 1)
    int nnz;          // Nombre de valeurs non nulles
    int nrows;        // Nombre de lignes de la matrice
} femIterativeSolver;

/* =========================================
   Conversion d'une matrice dense en CSR
   dense : matrice dense stockée en row-major
   nrows, ncols : dimensions de la matrice
   solver : structure où seront stockées les données CSR
   ========================================= */
void dense_to_csr(double *dense, int nrows, int ncols, femIterativeSolver *solver) {
    int i, j;
    // Comptage du nombre de valeurs non nulles
    int nnz = 0;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            if (fabs(dense[i * ncols + j]) > 1e-12)
                nnz++;
        }
    }
    solver->nnz = nnz;
    solver->nrows = nrows;

    // Allocation des tableaux CSR
    solver->values = (double*) malloc(nnz * sizeof(double));
    solver->col_index = (int*) malloc(nnz * sizeof(int));
    solver->row_ptr = (int*) malloc((nrows + 1) * sizeof(int));
    if (!solver->values || !solver->col_index || !solver->row_ptr) {
        fprintf(stderr, "Erreur d'allocation CSR\n");
        exit(EXIT_FAILURE);
    }

    int count = 0;
    solver->row_ptr[0] = 0;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            double val = dense[i * ncols + j];
            if (fabs(val) > 1e-12) {
                solver->values[count] = val;
                solver->col_index[count] = j;
                count++;
            }
        }
        solver->row_ptr[i+1] = count;
    }
}

/* =========================================
   Multiplication matrice-vecteur pour une matrice CSR
   y = A*x avec A stockée dans solver en CSR
   ========================================= */
void csr_matvec(const femIterativeSolver *solver, const double *x, double *y) {
    int i, j;
    for (i = 0; i < solver->nrows; i++) {
        y[i] = 0.0;
        for (j = solver->row_ptr[i]; j < solver->row_ptr[i+1]; j++) {
            y[i] += solver->values[j] * x[solver->col_index[j]];
        }
    }
}

/* =========================================
   Algorithme du Gradient Conjugué pour une matrice en format CSR
   Résout Ax = b pour une matrice SPD
   tol : tolérance du critère d'arrêt
   max_iter : nombre maximum d'itérations
   La solution est stockée dans solver->X
   ========================================= */
void conjugateGradientCSR(femIterativeSolver *solver, const double *b, double tol, int max_iter) {
    int n = solver->size;
    int i;
    
    // Allocation des vecteurs R, D et S si non déjà alloués
    solver->R = (double *) malloc(n * sizeof(double));
    solver->D = (double *) malloc(n * sizeof(double));
    solver->S = (double *) malloc(n * sizeof(double));
    solver->X = (double *) malloc(n * sizeof(double));
    if (!solver->R || !solver->D || !solver->S || !solver->X) {
        fprintf(stderr, "Erreur d'allocation des vecteurs du solveur\n");
        exit(EXIT_FAILURE);
    }

    // Initialisation : X = 0
    for (i = 0; i < n; i++) {
        solver->X[i] = 0.0;
    }

    // Calcul initial du résidu R = b - A*X, ici X = 0 donc R = b
    csr_matvec(solver, solver->X, solver->S); // S = A*X
    for (i = 0; i < n; i++) {
        solver->R[i] = b[i] - solver->S[i];
        solver->D[i] = solver->R[i];  // Initialisation de la direction
    }

    double rsold = 0.0;
    for (i = 0; i < n; i++) {
        rsold += solver->R[i] * solver->R[i];
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // S = A * D
        csr_matvec(solver, solver->D, solver->S);
        double pAp = 0.0;
        for (i = 0; i < n; i++) {
            pAp += solver->D[i] * solver->S[i];
        }
        double alpha = rsold / pAp;

        // Mise à jour de la solution X et du résidu R
        for (i = 0; i < n; i++) {
            solver->X[i] += alpha * solver->D[i];
            solver->R[i] -= alpha * solver->S[i];
        }
        double rsnew = 0.0;
        for (i = 0; i < n; i++) {
            rsnew += solver->R[i] * solver->R[i];
        }
        solver->error = sqrt(rsnew);
        solver->iter = iter + 1;

        // Critère d'arrêt
        if (solver->error < tol) {
            printf("Convergence atteinte en %d itérations, erreur = %e\n", iter+1, solver->error);
            return;
        }

        double beta = rsnew / rsold;
        for (i = 0; i < n; i++) {
            solver->D[i] = solver->R[i] + beta * solver->D[i];
        }
        rsold = rsnew;
    }
    printf("Nombre maximum d'itérations atteint, erreur = %e\n", solver->error);
}

/* =========================================
   Fonction main : création d'une matrice dense SPD,
   conversion en CSR, résolution par Gradient Conjugué et affichage de la solution.
   ========================================= */
int main(void) {
    // Exemple d'une matrice SPD 3x3
    int n = 3;
    double dense[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    // Vecteur second membre b
    double b[3] = { 1.0, 2.0, 3.0 };

    // Création et initialisation du solveur
    femIterativeSolver solver;
    solver.size = n;  // Taille du système
    // Conversion de la matrice dense en format CSR
    dense_to_csr(dense, n, n, &solver);

    // Paramètres de l'algorithme du Gradient Conjugué
    double tol = 1e-6;
    int max_iter = 1000;

    // Exécution de l'algorithme du Gradient Conjugué avec stockage CSR
    conjugateGradientCSR(&solver, b, tol, max_iter);

    // Affichage de la solution
    printf("Solution X =\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", solver.X[i]);
    }

    // Libération de la mémoire allouée pour le CSR et les vecteurs
    free(solver.values);
    free(solver.col_index);
    free(solver.row_ptr);
    free(solver.R);
    free(solver.D);
    free(solver.S);
    free(solver.X);

    return 0;
}
