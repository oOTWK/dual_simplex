/*** header file, dual simplex method for set-covering problems ***/

#ifndef DualSimplex_h
#define DualSimplex_h

int init_dual_simplex(int num_row, int num_col, int num_nonzero, int *row_sizes, int *col_sizes, int *costs, int **cols);
void destroy_simplex();
double solve_dual_simplex();

#endif /* DualSimplex_h */
