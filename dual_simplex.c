/*** 
Implementation of dual simplex method for set-covering problems (SCP)

Based on
 Koberstein, A. (2005). The dual simplex method, techniques for a fast and stable implementation. 
 Unpublished doctoral thesis, Universität Paderborn, Paderborn, Germany.

***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define BASIC 0
#define UPPER 1
#define LOWER 2

#define MAT(x,y) ((x) * xm + (y))

#define ZERO_TOL pow(10, -12)
#define D_FEA_TOL pow(10, -7)
#define P_FEA_TOL pow(10, -7)
#define R_FEA_TOL pow(10, -9)
#define DROP_TOL pow(10, -14)
#define PIVOT_TOL pow(10, -7)
#define TOL pow(10, -6)


typedef struct {
	int size;
	int *inds;
	double *values;
	char *is_in;
} ind_arr;

typedef struct {
	int ind;
	double value;
} pair;

typedef struct {
	int start;
	int ind;
} lbeg_elmt;

typedef struct {
	int ind;
	double value;
} l_node;

typedef struct u_node {
	union {int size; double value;};
	int rind;
	int cind;
	struct u_node *up;
	struct u_node *down;
	struct u_node *left;
	struct u_node *right;
} u_node;


static double obj_value;
static int xm, xn, xj;
static int num_fac;
static int *colwise_a;
static int *xjcp;
static double *xcost;
static int *l_bound;
static double *mat_b;
static double *xbxx, *xdjsc, *xbeta, *xbeta, *xpifs;
static int *xh, *xh_holder, *xh1, *xh2;
static char *xkey, *mark;

static ind_arr rho_r, tau, alpha_r, alpha_q, spike, xy;
static pair *al0, *au0;

static int lbeg_size, lbeg_size_f;
static l_node *eta_l;
static lbeg_elmt *lbeg;
static u_node *mat_u, *col_u_header, *row_u_header;

static l_node *eta_l_f;
static lbeg_elmt *lbeg_f;
static u_node *mat_u_f, *col_u_header_f, *row_u_header_f, *row_head_f;

static int *perm, *permback, *q_perm;
static int *row_entry, *row_link;


static void update_delta(int x_i);
static double compute_bfrt_theta(pair *al, pair *au, int al_size, int au_size);
static void add_u_node(u_node *u, int rind, int cind, double value, u_node *col_u_header, u_node *row_u_header);
static void remove_u_node(u_node *u, int rind, int cind, u_node *col_u_header, u_node *row_u_header);
static void f_tran(ind_arr *work, ind_arr *spike, int *last_p);
static void lu_fac();



#define MALLOC(var, type, size)     if ((var = (type) malloc(size)) == NULL) \
                                        { perror("Error malloc"); return -1; }
#define CALLOC(var, type, n, size)     if ((var = (type) calloc(n, size)) == NULL) \
                                        { perror("Error calloc"); return -1; }

int init_dual_simplex(int num_row, int num_col, int num_nonzero, int *row_sizes, int *col_sizes, int *costs, int **cols)
{
	int i, j;
	xm = num_row;
	xn = num_col;
	xj = xm + xn;

	MALLOC(colwise_a, int *, num_nonzero * sizeof(int))
	MALLOC(xjcp, int *, (xn+1) * sizeof(int))
	MALLOC(xcost, double *, xn * sizeof(double))
	MALLOC(l_bound, int *, xm * sizeof(int))
	MALLOC(xbxx, double *, xj * sizeof(double))
	MALLOC(xbeta, double *, xj * sizeof(double))
	MALLOC(xpifs, double *, xj * sizeof(double))

	CALLOC(xdjsc, double *, xj, sizeof(double))

	// index of basis
	MALLOC(xh1, int *, xm * sizeof(int))
	MALLOC(xh2, int *, xm * sizeof(int))
	xh = xh1;

	// status
	MALLOC(xkey, char *, xj * sizeof(char))

	CALLOC(mat_b, double *, xm * xm, sizeof(double))
	CALLOC(mark, char *, xm, sizeof(char))

	MALLOC(rho_r.inds, int *, xm * sizeof(int))
	CALLOC(rho_r.values, double *, xm, sizeof(double))
	rho_r.size = 0;
	tau = rho_r;

	MALLOC(alpha_r.inds, int *, xn * sizeof(int))
	MALLOC(alpha_r.values, double *, xn * sizeof(double))
	rho_r.size = 0;

	MALLOC(alpha_q.inds, int *, xm * sizeof(int))
	CALLOC(alpha_q.values, double *, xm, sizeof(double))
	alpha_q.size = 0;

	MALLOC(spike.inds, int *, xm * sizeof(int))
	CALLOC(spike.values, double *, xm, sizeof(double))
	spike.size = 0;

	MALLOC(xy.inds, int *, xm * sizeof(int))
	CALLOC(xy.values, double *, xm, sizeof(double))
	xy.size = 0;

	MALLOC(al0, pair *, xn * sizeof(pair))
	MALLOC(au0, pair *, xn * sizeof(pair))

	obj_value = 0;

	for (i = 0; i < num_row; i++) {
		l_bound[i] = - row_sizes[i];
	}

	for (i = 0, j = 0; i < num_col; i++) {
		xjcp[i] = j;
		xcost[i] = xdjsc[i] = ((double) costs[i]);
		memcpy(colwise_a + j, cols[i], col_sizes[i] * sizeof(int));
		j += col_sizes[i];
	}
	xjcp[i] = j;

	// init other vectors
	memset(xkey, LOWER, xn);
	for (i = 0; i < xm; i++) {
		xbxx[xn+i] = 0.0;
		xbeta[xn+i] = 1.0;
		xpifs[xn+i] = 1.0;
		xh[i] = xn + i;
		xkey[xn+i] = BASIC;
		mat_b[MAT(i,i)] = 1.0;
	}

	// init LU
	MALLOC(perm, int *, xm * sizeof(int))
	MALLOC(permback, int *, xm * sizeof(int))
	MALLOC(q_perm, int *, xm * sizeof(int))

	MALLOC(lbeg, lbeg_elmt *, xj * sizeof(lbeg_elmt))
	MALLOC(eta_l, l_node *, xm * xj * sizeof(l_node))
	lbeg_size = 0;
	lbeg[0].start = 0;


	CALLOC(mat_u, u_node *, xm * xm, sizeof(u_node))
	MALLOC(col_u_header, u_node *, xm * sizeof(u_node))
	MALLOC(row_u_header, u_node *, xm * sizeof(u_node))

	MALLOC(lbeg_f, lbeg_elmt *, xj * sizeof(lbeg_elmt))
	MALLOC(eta_l_f, l_node *, xm * xj * sizeof(l_node))
	lbeg_size_f = 0;
	lbeg_f[0].start = 0;


	CALLOC(mat_u_f, u_node *, xm * xm, sizeof(u_node))
	MALLOC(col_u_header_f, u_node *, xm * sizeof(u_node))
	MALLOC(row_head_f, u_node *, (xm+1) * sizeof(u_node))
	row_u_header_f = row_head_f + 1;
	row_head_f->rind = -1;

	MALLOC(row_entry, int *, xm * sizeof(int))
	MALLOC(row_link, int *, xm * sizeof(int))

	u_node *u;
	for (i = 0; i < xm; i++) {
		for (j = 0; j < xm; j++) {
			u = mat_u + MAT(i,j);
			u->rind = i;
			u->cind = j;
			u->up = NULL;

			u = mat_u_f + MAT(i,j);
			u->rind = i;
			u->cind = j;
			u->up = NULL;	
		}
	}

	for (i = 0; i < xm; i++) {
		u = col_u_header + i;
		u->up = u;
		u->down = u;
		u->rind = -1;
		u->cind = i;

		u = row_u_header + i;
		u->left = u;
		u->right = u; 
		u->rind = i;
		u->cind = -1;

		u = mat_u + MAT(i,i);
		add_u_node(u, i, i, 1.0, col_u_header, row_u_header);
		perm[i] = permback[i] = i;

		u = col_u_header_f + i;
		u->up = u;
		u->down = u;
		u->rind = -1;
		u->cind = i;

		u = row_u_header_f + i;
		u->left = u;
		u->right = u; 
		u->rind = i;
		u->cind = -1;
	}
	return 0;
}


void destroy_simplex()
{
	free(colwise_a);
	free(xjcp);
	free(xcost);
	free(l_bound);
	free(xbxx);
	free(xdjsc);
	free(xbeta);
	free(xpifs);

	free(xh1);
	free(xh2);
	free(xkey);
	free(mat_b);

	free(mark);
	free(rho_r.inds);
	free(rho_r.values);
	free(alpha_q.inds);
	free(alpha_q.values );
	free(spike.inds);
	free(spike.values);
	free(xy.inds);
	free(xy.values);

	free(al0);
	free(au0);
	free(perm);
	free(permback);
	free(q_perm);
	free(eta_l);
	free(mat_u);
	free(col_u_header);
	free(row_u_header);
	free(lbeg_f);
	free(eta_l_f);
	free(mat_u_f);
	free(col_u_header_f);
	free(row_head_f);
	free(row_entry);
	free(row_link);
}


double solve_dual_simplex()
{
	int i, j;
	int r_i, p_i, q_i = 0;
	double value, max_value;
	double alpha_rq = 0;
	double beta_r;
	double delta_zero, delta_one;
	int last;
	u_node *u;
	int refactor_freq;

	refactor_freq = 100 + xm / 200;
	if (refactor_freq > 2000) {
		refactor_freq = 2000;
	}

	int loop_c = 0;
	int fac_c = 0;
	num_fac = 0;

	while (1) {

		// pricing
		r_i = xm;
		max_value = -1;
		for (i = 0; i < xm; i++) {
			if (xpifs[xh[i]] != 0) { 
				value = xpifs[xh[i]] * xpifs[xh[i]] / xbeta[xh[i]];
				if (value > max_value) {
					max_value = value;
					r_i = i;
				}
			}
		}

		// optimal
		if (r_i == xm) {
			for (i = 0; i < xy.size; i++) {
				xy.values[xy.inds[i]] = 0.0;
			}

			xy.size = 0;
			for (i = 0; i < xm; i++) {
				if (xh[i] < xn) {
					xy.values[i] = xcost[xh[i]];
				}
			}

			for (i = 0; i < xm; i++) {
				j = perm[i];
				value = xy.values[j] / mat_u[MAT(j, j)].value;

				if (value >= ZERO_TOL || value <= - ZERO_TOL) {
					for (u = row_u_header[j].right; u->cind >= 0; u = u->right) {
						xy.values[u->cind] -= u->value * value;
					}

					xy.inds[xy.size] = j;
					xy.size++;
					xy.values[j] = value;
					mark[j] = 1;
				} else {
					xy.values[j] = 0.0;
				}		

			}
				
			for (i = lbeg_size-1; i >= 0; i--) {
				value = xy.values[lbeg[i].ind];
				if (value >= ZERO_TOL || value <= - ZERO_TOL) {
					for (j = lbeg[i].start; j < lbeg[i+1].start; j++) {
						if (mark[eta_l[j].ind] == 0) {
							mark[eta_l[j].ind] = 1;
							xy.inds[xy.size] = eta_l[j].ind;
							xy.size++;
						}
						xy.values[eta_l[j].ind] += eta_l[j].value * value;
					}
				} 
			}

			for (i = lbeg_size_f-1; i >= 0; i--) {
				value = xy.values[lbeg_f[i].ind];
				for (j = lbeg_f[i].start; j < lbeg_f[i+1].start; j++) {
					value += eta_l_f[j].value * xy.values[eta_l_f[j].ind];
				}

				if (value >= DROP_TOL || value <= - DROP_TOL) {
					if (mark[lbeg_f[i].ind] == 0) {
						mark[lbeg_f[i].ind] = 1;
						xy.inds[xy.size] = lbeg_f[i].ind;
						xy.size++;
					}
					xy.values[lbeg_f[i].ind] = value;	
				} else {
					xy.values[lbeg_f[i].ind] = 0.0;
				}
			}

			for (i = 0, j = 0; i < xy.size; i++) {
				value = xy.values[xy.inds[i]];
				mark[xy.inds[i]] = 0;

				if (value >= DROP_TOL || value <= - DROP_TOL) {
					xy.inds[j] = xy.inds[i];
					j++;
				} else {
					xy.values[xy.inds[i]] = 0.0;
				}
			}
			xy.size = j;

			printf("itr: %d\n", loop_c);
			return obj_value;
		} 

		fac_c++;
		if (fac_c > refactor_freq) {
			fac_c = 0;
			lu_fac();
			continue;
		}

		loop_c++;
		p_i = xh[r_i];
		delta_one = delta_zero = xpifs[p_i];


		// compute rho_r
		for (i = 0; i < rho_r.size; i++) {
			rho_r.values[rho_r.inds[i]] = 0.0;
		}

		rho_r.size = 0;
		rho_r.values[r_i] = 1.0;
		for (i = permback[r_i]; i < xm; i++) {
			j = perm[i];
			value = rho_r.values[j] / mat_u[MAT(j, j)].value;

			if (value >= ZERO_TOL || value <= - ZERO_TOL) {
				for (u = row_u_header[j].right; u->cind >= 0; u = u->right) {
					rho_r.values[u->cind] -= u->value * value;
				}

				rho_r.inds[rho_r.size] = j;
				rho_r.size++;
				rho_r.values[j] = value;
				mark[j] = 1;
			} else {
				rho_r.values[j] = 0.0;
			}
		}

		for (i = lbeg_size-1; i >= 0; i--) {
			value = rho_r.values[lbeg[i].ind];
			if (value >= ZERO_TOL || value <= - ZERO_TOL) {
				for (j = lbeg[i].start; j < lbeg[i+1].start; j++) {
					if (mark[eta_l[j].ind] == 0) {
						mark[eta_l[j].ind] = 1;
						rho_r.inds[rho_r.size] = eta_l[j].ind;
						rho_r.size++;
					}
					rho_r.values[eta_l[j].ind] += eta_l[j].value * value;
				}
			}
		}

		for (i = lbeg_size_f-1; i >= 0; i--) {
			value = rho_r.values[lbeg_f[i].ind];
			for (j = lbeg_f[i].start; j < lbeg_f[i+1].start; j++) {
				value += eta_l_f[j].value * rho_r.values[eta_l_f[j].ind];
			}

			if (value >= DROP_TOL || value <= - DROP_TOL) {
				if (mark[lbeg_f[i].ind] == 0) {
					mark[lbeg_f[i].ind] = 1;
					rho_r.inds[rho_r.size] = lbeg_f[i].ind;
					rho_r.size++;
				}
				rho_r.values[lbeg_f[i].ind] = value;	
			} else {
				rho_r.values[lbeg_f[i].ind] = 0.0;
			}
		}

		beta_r = 0;
		for (i = 0, j = 0; i < rho_r.size; i++) {
			value = rho_r.values[rho_r.inds[i]];
			mark[rho_r.inds[i]] = 0;
			if (value >= DROP_TOL || value <= - DROP_TOL) {
				rho_r.inds[j] = rho_r.inds[i];
				j++;
				beta_r += value * value;
			} else {
				rho_r.values[rho_r.inds[i]] = 0.0;
			}
		}
		rho_r.size = j;
		xbeta[p_i] = beta_r;


		// compute alpha_r (pivot row) this is packed array
		alpha_r.size = 0;
		int nonbasic_c = 0;
		for (i = 0; i < xn; i++) {
			value = 0;
			if (xkey[i] != BASIC) {
				nonbasic_c++;
				value = 0;
				for (j = xjcp[i]; j < xjcp[i+1]; j++) {
					value += rho_r.values[colwise_a[j]];
				}

				if (value >= DROP_TOL || value <= - DROP_TOL) {
					alpha_r.inds[alpha_r.size] = i;
					alpha_r.values[alpha_r.size] = value;
					alpha_r.size++;
				}				
			}
		}
		for (i = xn; i < xj && nonbasic_c < xn; i++) {	
			if (xkey[i] != BASIC) {
				nonbasic_c++;
				value = rho_r.values[i-xn];

				if (value >= DROP_TOL || value <= -DROP_TOL) {
					alpha_r.inds[alpha_r.size] = i;
					alpha_r.values[alpha_r.size] = value;
					alpha_r.size++;
				}	
			}
		}


		// ========================== ratio test ==========================
		double theta_D;
		double theta;

		int j1, j2;
		j1 = j2 = 0;
		if (xpifs[p_i] > 0) {
			for (i = 0; i < alpha_r.size; i++) {
				if (alpha_r.values[i] > PIVOT_TOL && xkey[alpha_r.inds[i]] == LOWER) {
					al0[j1].value = alpha_r.values[i];
					al0[j1].ind = alpha_r.inds[i];
					j1++;
				} else if (alpha_r.values[i] < - PIVOT_TOL && xkey[alpha_r.inds[i]] == UPPER) {
					au0[j2].value = alpha_r.values[i];
					au0[j2].ind = alpha_r.inds[i];
					j2++;
				}
			}

		} else if (xpifs[p_i] < 0){
			for (i = 0; i < alpha_r.size; i++) {
				if (alpha_r.values[i] < - PIVOT_TOL && xkey[alpha_r.inds[i]] == LOWER) {
					al0[j1].value = - alpha_r.values[i];
					al0[j1].ind = alpha_r.inds[i];
					j1++;
				} else if (alpha_r.values[i] > PIVOT_TOL && xkey[alpha_r.inds[i]] == UPPER) {
					au0[j2].value = - alpha_r.values[i];
					au0[j2].ind = alpha_r.inds[i];
					j2++;
				}
			}
		}

		theta = compute_bfrt_theta(al0, au0, j1, j2);

		double max_value = -1.0;
		for (i = 0; i < j1; i++) {
			if (xdjsc[al0[i].ind] / al0[i].value <= theta) {
				value = al0[i].value;
				if (value > max_value) {
					max_value = value;
					q_i = al0[i].ind;
					alpha_rq = al0[i].value;
				}	
			}
		}

		for (i = 0; i < j2; i++) {
			if (xdjsc[au0[i].ind] / au0[i].value <= theta) {
				value = - au0[i].value;
				if (value > max_value) {
					max_value = value;
					q_i = au0[i].ind;
					alpha_rq = au0[i].value;
				}
			}
		}

		if (xpifs[p_i] < 0) {
			alpha_rq *= -1.0;
		}

		theta_D = xdjsc[q_i] / alpha_rq;

		if (xpifs[p_i] > 0) {
			if (theta_D < 0) {
				theta_D = 0.0;
			}
		} else {
			if (theta_D > 0) {
				theta_D = 0.0;
			}
		}
		// ============================ end of ratio test ==========================


		// compute alpha_q (pivot column)
		for (i = 0; i < alpha_q.size; i++) {
			alpha_q.values[alpha_q.inds[i]] = 0.0;
		}

		alpha_q.size = 0;
		if (q_i < xn) {
			for (i = xjcp[q_i]; i < xjcp[q_i+1]; i++) {
				alpha_q.values[colwise_a[i]] = 1.0;
				alpha_q.inds[alpha_q.size] = colwise_a[i];
				alpha_q.size++;
				mark[colwise_a[i]] = 1;
			}
		} else {
			alpha_q.values[q_i - xn] = 1.0;
			alpha_q.inds[0] = q_i - xn;
			alpha_q.size = 1;
			mark[q_i - xn] = 1;
		}

		for (i = 0; i < spike.size; i++) {
			spike.values[spike.inds[i]] = 0.0;
		}

		spike.size = 0;
		f_tran(&alpha_q, &spike, &last);

		// check alpha_rq diff
		double threshold;
		threshold = alpha_q.values[r_i]; if (threshold < 0) { threshold *= -1.0; }
		threshold = pow(10, -9) * (1 + threshold);
		if (alpha_q.values[r_i] - alpha_rq > threshold || alpha_q.values[r_i] - alpha_rq < -threshold) {
			lu_fac();
			continue;
		}


		/********* update according to p <-> q *********/
		// compute tau
		tau.size = rho_r.size;
		for (i = 0; i < tau.size; i++) {
			mark[tau.inds[i]] = 1;
		}
		
		f_tran(&tau, NULL, NULL);
		rho_r.size = tau.size;


		// update xbeta
		beta_r = beta_r / pow(alpha_q.values[r_i], 2);
		double k = - 2.0 / alpha_q.values[r_i];
		for (i = 0; i < alpha_q.size; i++) {
			j = alpha_q.inds[i];
			xbeta[xh[j]] = xbeta[xh[j]] + alpha_q.values[j] * (alpha_q.values[j] * beta_r + k * tau.values[j]);

			if (xbeta[xh[j]] < 0.0001) {
				xbeta[xh[j]] = 0.00001;
			}
		}
		xbeta[p_i] = 0.0;
		xbeta[q_i] = beta_r;


		// update xdjsc (reduced cost)
		for (i = 0; i < alpha_r.size; i++) {
			j = alpha_r.inds[i];
			xdjsc[j] -= theta_D * alpha_r.values[i];
		}

		xdjsc[p_i] = - theta_D;
		xdjsc[q_i] = 0.0;

		double theta_P;
		theta_P = delta_zero / alpha_q.values[r_i];

		for (i = 0; i < alpha_q.size; i++) {
			j = xh[alpha_q.inds[i]];
			xbxx[j] -= theta_P * alpha_q.values[alpha_q.inds[i]];
			update_delta(j);
		}
		if (q_i >= xn) {
			if (xkey[q_i] == UPPER) {
				xbxx[q_i] = -1.0;
			} else if (xkey[q_i] == LOWER) {
				xbxx[q_i] = l_bound[q_i - xn];
			} else {
				printf("ERROR xkey in updating x_q\n");
			}
		} else {
			if (xkey[q_i] == UPPER) {
				xbxx[q_i] = 1.0;
			} else if (xkey[q_i] == LOWER) {
				xbxx[q_i] = 0.0;
			} else {
				printf("ERROR xkey in updating x_q\n");
			}
		}
		xbxx[q_i] += theta_P;
		update_delta(q_i);

		obj_value += theta_D * delta_zero; 


		// update xkey, xh, b
		xkey[p_i] = (delta_one > 0 ? UPPER : LOWER);

		if ((xdjsc[p_i] < 0 && xkey[p_i] == LOWER) || (xdjsc[p_i] > 0 && xkey[p_i] == UPPER)) {
			printf("wrong xkey %d, d %f, delta_0 %f\n", xkey[p_i], xdjsc[p_i], delta_zero);
			exit(0);
		}

		xkey[q_i] = BASIC;
		xh[r_i] = q_i;


		// compute permuted spike
		u_node  *v, *u2;
		double multiplier;
		double u_rr =  mat_u[MAT(r_i, r_i)].value;

		// remove r-th col from u
		for (u = col_u_header[r_i].down; u->rind >= 0; u = u->down) {
			remove_u_node(u, u->rind, u->cind, col_u_header, row_u_header);
		}

		// add spike into u at position r
		for (i = 0; i < spike.size; i++) {
			value = spike.values[spike.inds[i]];

			add_u_node(mat_u + MAT(spike.inds[i], r_i),
				spike.inds[i], r_i,
				spike.values[spike.inds[i]],
				col_u_header, row_u_header);
		}


		// update U and create L
		int eta_l_ind;
		eta_l_ind = lbeg[lbeg_size].start;

		for (j = permback[r_i] + 1; j <= last; j++) {
			i = perm[j];
			u = mat_u + MAT(r_i, i);

			if (u->value >= ZERO_TOL || u->value <= - ZERO_TOL) {
				multiplier = - u->value / mat_u[MAT(i,i)].value;

				eta_l[eta_l_ind].value = multiplier;
				eta_l[eta_l_ind].ind = i;
				eta_l_ind++;

				for (v = row_u_header[i].right; v->cind >= 0; v = v->right) {
					u2 = mat_u + MAT(r_i, v->cind);
					value = u2->value + multiplier * v->value;
					if (value >= DROP_TOL || value <= - DROP_TOL) {
						if (u2->up == NULL) {
							add_u_node(u2, r_i, v->cind, value, col_u_header, row_u_header);
						} else {
							u2->value = value;
						}
					} else {
						u2->value = 0.0;
					}
				}				
			}
			u->value = 0.0;
		}

		if (eta_l_ind > lbeg[lbeg_size].start) {
			lbeg[lbeg_size].ind = r_i;
			lbeg_size++;
			lbeg[lbeg_size].start = eta_l_ind;
		}

		// remove zero node from r-rh row of u 
		for (u = row_u_header[r_i].right; u->cind >= 0; u = u->right) {
			if (u->value < DROP_TOL && u->value > - DROP_TOL) {
				remove_u_node(u, u->rind, u->cind, col_u_header, row_u_header);
			}
 		}  

		// update permutation
 		for (j = permback[r_i] + 1; j <= last; j++) {
 			i = perm[j];
 			permback[i] = j - 1;
 			perm[j - 1] = i;
 		}
 		permback[r_i] = last;
 		perm[last] = r_i;


 		// LU update stability test
 		double tol;
 		tol = alpha_q.values[r_i];
 		if (tol < 0) {
 			tol *= -1.0;
 		}
 		tol *= pow(10, -6);

 		value = mat_u[MAT(perm[last], perm[last])].value / u_rr;

 		// TODO
 		if (alpha_q.values[r_i] - value > tol || alpha_q.values[r_i] - value < - tol) {
 			printf("ERROR basis change should be aborted %f %f\n", alpha_q.values[r_i], value);
 			//exit(0);
 		}
	}

	return obj_value;
}


static void update_delta(int x_i)
{
	double x_value = xbxx[x_i];
	double *delta_p = xpifs + x_i;

	if (x_i < xn) {
		if (x_value < - P_FEA_TOL) {
			*delta_p = x_value;
		} else if (x_value > 1.0 + P_FEA_TOL + R_FEA_TOL) {
			*delta_p = x_value - 1.0;
		} else { // primal feasible
			*delta_p = 0.0;
		}
	} else {
		if (x_value > -1.0 + R_FEA_TOL + P_FEA_TOL) {
			*delta_p = x_value + 1.0;
		} else if (x_value < l_bound[x_i - xn] + l_bound[x_i - xn]*R_FEA_TOL - P_FEA_TOL) {

			*delta_p = x_value - l_bound[x_i - xn];
		}else {
			*delta_p = 0.0;
		}
	}
}


static double compute_bfrt_theta(pair *al, pair *au, int al_size, int au_size)
{
	double theta;
	int i;
	double value;

	theta = -1.0;
	for (i = 0; i < al_size; i++) {
		value = (xdjsc[al[i].ind] + D_FEA_TOL) / al[i].value;
		if (theta == -1.0 || value < theta) {
			theta = value;
		}
	}
	for (i = 0; i < au_size; i++) {
		value = (xdjsc[au[i].ind] - D_FEA_TOL) / au[i].value;
		if (theta == -1.0 || value < theta) {
			theta = value;
		}
	}

	// TODO: should deal with this
	if (theta < 0) {
		printf("WARNING negative theta %f\n", theta);
		theta = 0;
		//exit(0);
	}

	return theta;
}


static void f_tran(ind_arr *work, ind_arr *spike, int *last_p)
{
	int i, j;
	double value;
	int last;
	u_node *u;

	// FTranL-F
	for (i = 0; i < lbeg_size_f; i++) {
		value = work->values[lbeg_f[i].ind];
		if (value >= ZERO_TOL || value <= - ZERO_TOL) {
			for (j = lbeg_f[i].start; j < lbeg_f[i+1].start; j++) {
				if (mark[eta_l_f[j].ind] == 0) {
					mark[eta_l_f[j].ind] = 1;
					work->inds[work->size] = eta_l_f[j].ind;
					work->size++;
				}
				work->values[eta_l_f[j].ind] += eta_l_f[j].value * value;
			}
		}
	}

	// FTranL-U
	for (i = 0; i < lbeg_size; i++) {
		value = work->values[lbeg[i].ind];
		for (j = lbeg[i].start; j < lbeg[i+1].start; j++) {
			value += eta_l[j].value * work->values[eta_l[j].ind];
		}
		if (value >= DROP_TOL || value <= - DROP_TOL) {
			if (mark[lbeg[i].ind] == 0) {
				mark[lbeg[i].ind] = 1;
				work->inds[work->size] = lbeg[i].ind;
				work->size++;
			}
			work->values[lbeg[i].ind] = value;	
		} else {
			work->values[lbeg[i].ind] = 0.0;
		}
	}

	last = 0;
	for (i = 0; i < work->size; i++) {
		value = work->values[work->inds[i]];
		mark[work->inds[i]] = 0;
		if (value >= DROP_TOL || value <= - DROP_TOL) {
			if (permback[work->inds[i]] > last) {
				last = permback[work->inds[i]];
			}

			if (spike != NULL) {
				spike->values[work->inds[i]] = value;
				spike->inds[spike->size] = work->inds[i];
				spike->size++;
			}
		} 
		else {
			work->values[work->inds[i]] = 0.0;
		}
	}

	if (last_p != NULL) {
		*last_p = last;
	}

	work->size = 0;
	for (i = last; i >= 0; i--) {
		j = perm[i];
		value = work->values[j] / mat_u[MAT(j, j)].value;
		if (value >= ZERO_TOL || value <= - ZERO_TOL) {
			for (u = col_u_header[j].down; u->rind >= 0; u = u->down) {
				work->values[u->rind] -= u->value * value;
			}
			work->inds[work->size] = j;
			work->size++;
			work->values[j] = value;
		}
		 else {
			work->values[j] = 0.0;
		}
	}
}


static void add_u_node(u_node *u, int rind, int cind, double value, 
					   u_node *col_u_header, u_node *row_u_header)
{
	u_node *col_head = col_u_header + cind;
	u_node *row_head = row_u_header + rind;

	u->up = col_head->up;
	u->down = col_head;
	col_head->up->down = u;
	col_head->up = u;

	u->left = row_head->left;
	u->right = row_head;
	row_head->left->right = u;
	row_head->left = u;

	u->value = value;
}


static void remove_u_node(u_node *u, int rind, int cind, 
						  u_node *col_u_header, u_node *row_u_header)
{
	u->down->up = u->up;
	u->up->down = u->down;

	u->right->left = u->left;
	u->left->right = u->right;

	u->value = 0.0;
	u->up = NULL;
}


static void lu_fac()
{
	num_fac++;
	int i, j, k;

	u_node *u, *r, *v;
	u_node *pivot_u;

	// clean up and init LU fac
	lbeg_size_f = 0;
	lbeg_f[0].start = 0;
	for (i = 0; i < xm; i++) {
		for (u = row_u_header_f[i].right; u->cind >= 0; u = u->right) {
			u->up = NULL;
			u->value = 0.0;
		}
	}

	for (i = 0; i < xm; i++) {
		u = col_u_header_f + i;
		u->up = u;
		u->down = u;
		u->size = -1;

		u = row_u_header_f + i;
		u->left = u;
		u->right = u; 
		u->size = -1;

		(u-1)->down = u;
		u->up = u-1;
	}

	row_head_f->up = row_head_f+xm;
	(row_head_f+xm)->down = row_head_f;


	for (i = 0; i < xm; i++) {
		j = xh[i];

		if (j < xn) {
			for (k = xjcp[j]; k < xjcp[j+1]; k++) {
				u = mat_u_f + MAT(colwise_a[k], i);
				add_u_node(u, colwise_a[k], i, 1.0, col_u_header_f, row_u_header_f);
				row_u_header_f[colwise_a[k]].size++;
				col_u_header_f[i].size++;
			}
		} else {
			u = mat_u_f + MAT(j-xn, i);
			add_u_node(u, j-xn, i, 1.0, col_u_header_f, row_u_header_f);
			row_u_header_f[j-xn].size++;
			col_u_header_f[i].size++;			
		}
	}

	// init row_entry, row_link
	for (i = 0; i < xm; i++) {
		row_entry[i] = -1;
	}
	for (r = row_head_f->up; r != row_head_f; r = r->up) {
		j = r->size;
		if (row_entry[j] < 0) {
			row_link[r->rind] = -1;
		} else {
			row_link[r->rind] = row_entry[j];
		}
		row_entry[j] = r->rind;
	}

	double max_value, value;
	int min_score, score_threshold, score, min_local_score;
	int *prev_row_link = row_entry;

	for (k = 0; k < xm ; k++) {

		// find pivot
		if ((j = row_entry[0]) >= 0) {
			pivot_u = row_u_header_f[j].right;
			row_entry[0] = row_link[j];
		} else {
			min_score = xm * xm;
			pivot_u = NULL;
			for (i = 0, j = -1; ; prev_row_link = row_link + j, j = *prev_row_link) {
				while (j < 0 && i < xm) {
					i++;
					prev_row_link = row_entry + i;
					j = *prev_row_link;
					
				}
				if (i == xm) {
					break;
				}

				// find pivot in j-th row
				v = NULL;
				max_value = 0.0;
				score_threshold = (min_score / i) + 1;
				if (score_threshold * i == min_score) {
					score_threshold--;
				}
				min_local_score = score_threshold;

				for (u = row_u_header_f[j].right; u->cind >= 0; u = u->right) {
					if ((score = col_u_header_f[u->cind].size) == 0) {
						pivot_u = u;
						*prev_row_link = row_link[j];
						i = xm;
						break;
					}

					if (score < min_local_score) {
						min_local_score = score;
						v = u; // candidate pivot_u
					}

					value = (u->value < 0 ? - u->value : u->value);
					if (value > max_value) {
						max_value = value;
					}
				}

				if (i == xm) {
					break;
				}

				if (v != NULL) {
					if ((v->value < 0 ? - v->value : v->value) >= 0.01 * max_value) {
						pivot_u = v;
						min_score = i * min_local_score;
					} else {
						for (u = row_u_header_f[j].right; u->cind >= 0; u = u->right) {

							if ((score = col_u_header_f[u->cind].size) < score_threshold &&
							 (u->value < 0 ? - u->value : u->value) >= 0.01 * max_value) {

								score_threshold = score;
								pivot_u = u;
							}
						}

						if (pivot_u->rind == j) {
							min_score = i * score_threshold;
						}
					}
				}
			}
		}

		// eliminate
		int pivot_c, pivot_r;
		double multiplier;
		int eta_l_ind;

		pivot_c = pivot_u->cind;
		pivot_r = pivot_u->rind;

		// permutation
		permback[pivot_r] = k;

		perm[k] = pivot_r;
		q_perm[pivot_c] = k;

		// remove pivot row
		u = row_u_header_f + pivot_r;
		u->down->up = u->up;
		u->up->down = u->down;
		u->up = NULL; // necessary?

		// alter col_size of pivot row
		for (u = row_u_header_f[pivot_r].right; u->cind >= 0; u = u->right) {
			col_u_header_f[u->cind].size--;
		}

		if (col_u_header_f[pivot_c].size == -1) {
			continue;
		}

		// init row_entry
		for (i = 0; i < xm; i++) {
			row_entry[i] = -1;
		}

		eta_l_ind = lbeg_f[lbeg_size_f].start;
		for (r = row_head_f->up; r != row_head_f; r = r->up) {
			u = mat_u_f + MAT(r->rind, pivot_c);
			if (u->up != NULL) {
				value = u->value;
				remove_u_node(u, u->rind, u->cind, col_u_header_f, row_u_header_f);
				if (value >= ZERO_TOL || value <= - ZERO_TOL) {
					multiplier = - value / pivot_u->value;

					eta_l_f[eta_l_ind].value = multiplier;
					eta_l_f[eta_l_ind].ind = u->rind;
					eta_l_ind++;

					for (u = row_u_header_f[pivot_r].right; u->cind >= 0; u = u->right) {
						if (u->cind == pivot_c) {
							continue;
						}

						v = mat_u_f + MAT(r->rind, u->cind);
						value = v->value + multiplier * u->value;

						if (value < DROP_TOL && value > - DROP_TOL) {
							if (v->up != NULL) {
								remove_u_node(v, v->rind, v->cind, col_u_header_f, row_u_header_f);
								row_u_header_f[v->rind].size--;
								col_u_header_f[v->cind].size--;
							} else if (v->value != 0.0) {
								printf("ERROR non zero u value %f at %d\n", u->value, k); exit(0);
							}
						} else {
							if (v->up == NULL) {
								add_u_node(v, v->rind, v->cind, value, col_u_header_f, row_u_header_f);
								row_u_header_f[v->rind].size++;
								col_u_header_f[v->cind].size++;		
							} else {
								v->value = value;
							}		
						}
					}
				}
			}

			// insert 
			j = r->size;
			if (row_entry[j] < 0) {
				row_link[r->rind] = -1;
			} else {
				row_link[r->rind] = row_entry[j];
			}
			row_entry[j] = r->rind;			
		}

		if (eta_l_ind > lbeg_f[lbeg_size_f].start) {
			lbeg_f[lbeg_size_f].ind = pivot_r;
			lbeg_size_f++;
			lbeg_f[lbeg_size_f].start = eta_l_ind;
		}
	}

	if (xh == xh1) {
		xh_holder = xh1;
		xh = xh2;
	} else {
		xh_holder = xh2;
		xh = xh1;
	}

	lbeg_size = 0;
	lbeg[0].start = 0;
	for (i = 0; i < xm; i++) {
		for (u = row_u_header[i].right; u->cind >= 0; u = u->right) {
			u->up = NULL;
			u->value = 0.0;
		}
	}

	for (i = 0; i < xm; i++) {
		u = col_u_header + i;
		u->up = u;
		u->down = u;

		u = row_u_header + i;
		u->left = u;
		u->right = u; 
	}

	for (i = 0; i < xm; i++) {
		xh[perm[q_perm[i]]] = xh_holder[i];
		for (u = col_u_header_f[i].down; u->rind >= 0; u = u->down) {
			v = mat_u + MAT(u->rind, perm[q_perm[i]]);
			add_u_node(v, v->rind, v->cind, u->value, col_u_header, row_u_header);
		}
	}

	for (i = 0; i < xy.size; i++) {
		xy.values[xy.inds[i]] = 0.0;
	}

	xy.size = 0;

	for (i = 0; i < xm; i++) {
		if (xh[i] < xn) {
			xy.values[i] = xcost[xh[i]];
		}
	}

	for (i = 0; i < xm; i++) {
		j = perm[i];
		value = xy.values[j] / mat_u[MAT(j, j)].value;
		if (value >= ZERO_TOL || value <= - ZERO_TOL) {
			for (u = row_u_header[j].right; u->cind >= 0; u = u->right) {
				xy.values[u->cind] -= u->value * value;
			}
			xy.inds[xy.size] = j;
			xy.size++;
			xy.values[j] = value;
			mark[j] = 1;
		} else {
			xy.values[j] = 0.0;
		}		

	}

	for (i = lbeg_size_f-1; i >= 0; i--) {
		value = xy.values[lbeg_f[i].ind];
		for (j = lbeg_f[i].start; j < lbeg_f[i+1].start; j++) {
			value += eta_l_f[j].value * xy.values[eta_l_f[j].ind];
		}

		if (value >= DROP_TOL || value <= - DROP_TOL) {
			if (mark[lbeg_f[i].ind] == 0) {
				mark[lbeg_f[i].ind] = 1;
				xy.inds[xy.size] = lbeg_f[i].ind;
				xy.size++;
			}
			xy.values[lbeg_f[i].ind] = value;	
		} else {
			xy.values[lbeg_f[i].ind] = 0.0;
		}
	}

	for (i = 0, j = 0; i < xy.size; i++) {
		value = xy.values[xy.inds[i]];
		mark[xy.inds[i]] = 0;

		if (value >= DROP_TOL || value <= - DROP_TOL) {
			xy.inds[j] = xy.inds[i];
			j++;
		} else {
			xy.values[xy.inds[i]] = 0.0;
		}
	}
	xy.size = j;

	for (i = 0; i < xj; i++) {
		if (i < xn) {
			value = 0.0;
			for (j = xjcp[i]; j < xjcp[i+1]; j++) {
				value += xy.values[colwise_a[j]];
			}
			value = xcost[i] - value;			
		} else {
			value = - xy.values[i - xn];	
		}

		// recompute
		xdjsc[i] = value;
		if (xkey[i] != BASIC) {
			if (value >= D_FEA_TOL) {
				xkey[i] = LOWER;
			} else if (value <= - D_FEA_TOL) {
				xkey[i] = UPPER;
			}
		}
	}

	for (i = 0; i < xy.size; i++) {
		xy.values[xy.inds[i]] = 0.0;
	}

	xy.size = 0;

	for (i = 0; i < xj; i++) {
		if (xkey[i] == BASIC) 
			continue;

		if (i < xn) {
			if (xkey[i] == UPPER) {
				for (j = xjcp[i]; j < xjcp[i+1]; j++) {
					xy.values[colwise_a[j]] -= 1.0;
				}
			}
		} else {
			if (xkey[i] == UPPER) {
				xy.values[i - xn] += 1.0;
			} else {
				xy.values[i - xn] -= l_bound[i-xn];
			}
		}
	}

	for (i = 0; i < xm; i++) {
		if (xy.values[i] != 0.0) {
			mark[i] = 1;
			xy.inds[xy.size] = i;
			xy.size++;
		}
	}

	f_tran(&xy, NULL, NULL);

	for (i = 0; i < xm; i++) {
		// recompute
		xbxx[xh[i]] = xy.values[i];
		update_delta(xh[i]);
	}
}
