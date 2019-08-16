/*
PageRank_iterations should implement the iterative procedure of the PageRank algorithm.
The values of the damping constant d and the convegence threshold epsilon should be
among the input arguments. The converged PageRank score vector x should be among the
output arguments
*/

void spmm(int nodes, double *val, int *col_idx, int *row_ptr,
          double *x_old, double *spmm_result, double d){

  for (int i = 0; i < nodes; i++){
    spmm_result[i] = 0;
    for (j = row_ptr[i]; j < row_ptr[i+1] - 1; j++){
        spmm_result[i] += val[j] * x_old[col_idx[j]];
    spmm_result[i] *= d;
    } //end of inner loop

  }   // end of outer loop

}     // end of function

void dangling_PageRank(double *x_old, int *dangling_idx){
  int n = sizeof(dangling_idx)/sizeof(int)
  int sum = 0;

  for (int i = 0; i < n; i++) {
    sum += x_old[dangling_idx[i]];
  }
  return sum;
}

void PageRank_iterations(double *val, int *col_idx,
                         int *row_ptr, int *dangling_idx,
                         int nodes, double d) {
  int i;
  double W_sum;
  double *temp_ptr;


  double *x_old = (double*)malloc(nodes*sizeof(double));
  double *x_new = (double*)malloc(nodes*sizeof(double));
  double *spmm_result = (double*)malloc(nodes*sizeof(double));

  double one_over_N = 1.0/nodes;
  for (i = 0; i < nodes; i++){
    x_old[i] = one_over_N;
    x_new[i] = 1;
  }

  for (i = 0; i < nodes; i++){
    W_sum = dangling_PageRank(x_old, dangling_idx);
    spmm(nodes, val, col_idx, row_ptr, x_old, spmm_result, d);
    x_new[i] = ((1 - d + d*W_sum)*one_over_N) + spmm_result[i];
  }
  //Pointer swapping
  temp_ptr = x_new;
  x_new = x_old;
  x_old = temp_ptr;

  for (i = 0; i < nodes; i++) printf("%d  |  %lf\n", i, x_old[i]);

}

int main(){
  //argv[1] - firt argument


}
