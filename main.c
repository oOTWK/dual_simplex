/*** main.c for dual simplex method test on set-covering problems ***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dual_simplex.h"


#define FILE_FORMAT_ERR    fprintf(stderr, "Error: wrong SCP file format\n")
    
#define GETLINE(buf, buf_size, fp)  if (getline(&buf, &buf_size, fp) == -1) \
                                        { FILE_FORMAT_ERR; return -1; }

#define STR_TOKEN(token, s, delim)  if ((token = strtok(s, delim)) == NULL) \
                                        { FILE_FORMAT_ERR; return -1; }

#define MALLOC(var, type, size)     if ((var = (type) malloc(size)) == NULL) \
                                        { perror("Error malloc"); return -1; }



int read_input(char *filename)
{
    int i, j;
    int num_row, num_col, num_nonzero;
    int *costs;
	int *col_sizes;
	int *row_sizes;

    FILE *fp;
    char *buf = NULL;
    size_t buf_size = 0;
    char *token;

    if ((fp = fopen(filename, "r")) == NULL) { 
        perror("Error opening file"); return -1; 
    }

    // read first line: the number of row and the number of col
    GETLINE(buf, buf_size, fp)
    STR_TOKEN(token, buf, " ")
    num_row = atoi(token);
    STR_TOKEN(token, NULL, " ")
    num_col = atoi(token);

    if ((costs = (int *) malloc(num_col * sizeof(int))) == NULL) { 
        perror("Error malloc"); return -1;
    }

    // read cost vector
    for (i = 0, token = NULL; i < num_col; token = strtok(NULL, " ")) {
        if (token == NULL) {
            GETLINE(buf, buf_size, fp)
            STR_TOKEN(token, buf, " ")
        }

        if (*token != '\n') {
            costs[i++] = atoi(token);
        }
    }

    // temporary constraint matrix
    int **cols;
    MALLOC(cols, int **, num_col * sizeof(int *))
    for (i = 0; i < num_col; i++) {
        MALLOC(cols[i], int *, num_row * sizeof(int))
    }

    MALLOC(col_sizes, int *, num_col * sizeof(int))
    memset(col_sizes, 0, num_col * sizeof(int));
    MALLOC(row_sizes, int *, num_row * sizeof(int))
    memset(row_sizes, 0, num_row * sizeof(int));

    // read rows (constraints)
    int col_idx;
    num_nonzero = 0;

    for (i = 0; i < num_row; i++) {
        GETLINE(buf, buf_size, fp)
        STR_TOKEN(token, buf, " ")
        row_sizes[i] = atoi(token);
        num_nonzero += row_sizes[i];

        for (j = 0, token = NULL; j < row_sizes[i]; token = strtok(NULL, " ")) {
            if (token == NULL) {
                GETLINE(buf, buf_size, fp)
                STR_TOKEN(token, buf, " ")
            }

            if (*token != '\n') {
                if ((col_idx = atoi(token) - 1) < 0) {
                    FILE_FORMAT_ERR; return -1;
                }
                cols[col_idx][col_sizes[col_idx]] = i;
                col_sizes[col_idx]++;
                j++;
            }
        }
    }

    if (init_dual_simplex(num_row, num_col, num_nonzero, row_sizes, col_sizes, costs, cols)) 
    	return -1;

    // free temporary constraint matrix
    for (i = 0; i < num_col; i++) {
        free(cols[i]);
    }

    return 0;
}


int main(int argc, char *argv[])
{
	clock_t begin_t, end_t;
	double LP_soln;

	if (argc != 2) {
		fprintf(stderr, "usage: %s input_file\n", argv[0]);
		exit(1);
	}

	begin_t = clock();
	if (read_input(argv[1])) return 1;
	LP_soln = solve_dual_simplex();
	destroy_simplex();
	end_t = clock();

	printf("LP soln: %f\n", LP_soln);
	printf("CPU time %.3fs\n", (double) (end_t - begin_t) / CLOCKS_PER_SEC);
	
	return 0;
}
