
void read_graph_file(FILE *infile, double **val_ptr, int **col_idx_ptr,
              int **row_ptr_ptr, int **dangling_idx_ptr, int nodes, int edges) {

  int temp1, temp2;
  int i, j, k;

  char buffer[100]; // Array for storing lines that we have to skip

  int *FromNode = (int*)malloc(edges*sizeof(int)); //array for storing the FromNodeId-values
  int *ToNode = (int*)malloc(edges*sizeof(int)); // array for storing the ToNodeId-values
  int *num_out_links = (int*)malloc(nodes*sizeof(int)); // no. of outbound links from each webpage

  // Collecting the FromeNodeId and ToNodeId data
  i = 0;
  while (fgets(buffer, sizeof(buffer), infile) != NULL){
    sscanf(buffer, "%d %d", &temp1, &temp2);
    FromNode[i] = temp1;
    ToNode[i] = temp2;
    i++;
  } // end of while-loop
  fclose(infile);


  //Collecting total number of FromNode links for each website, and total number of nonzero values
  int nz = 0, z = 0; // no. of nonzero elements and number of zero elements
  /*The z-counter is for sizing the danging-webpages index array later*/
  int temp_counter; // Counter k


  #pragma omp parallel for private(i,j) reduction(+:nz)
  for (i = 0; i < nodes; i++){
    num_out_links[i] = 0; // Initialize L;

    for (j = 0; j < edges; j++){
      if ((FromNode[j] == i) && (ToNode[j] != FromNode[j])) {
      num_out_links[i] += 1;
      nz++;
      } // end of if-statement
    } // end of inner loop
  if (num_out_links[i] == 0) z++;
  } // end of outer loop

  if (z != 0){ // If there are any dangling webpages
    dangling_idx = (int*) realloc(dangling_idx, z*sizeof(int));
    temp_counter = 0;
    for (i = 0; i < nodes; i++){
      if (num_out_links[i] == 0){
        dangling_idx[temp_counter] = i;
        temp_counter++;
      } // end of inner if-statement
    }   // end of inner loop
  }     // end of outer if-statement
  else free(dangling_idx);


  /* Building the hyperlink matrix in CRS-format */
  *val_ptr     = (double*) malloc(nz*sizeof(double));
  *col_idx_ptr = (int*)    malloc(nz*sizeof(int));
  *row_ptr_ptr = (int*)    malloc((nodes+1)*sizeof(int));


  temp_counter = 0;
  for (i = 0; i < nodes; i++){
    row_ptr[i] = temp_counter;
    //printf("%d\n", row_ptr[i]);
    for (j = 0; j < edges; j++){

      if ((ToNode[j] == i) && (ToNode[j] != FromNode[j])){
        val[temp_counter] = 1.0/num_out_links[FromNode[j]];
        col_idx[temp_counter] = FromNode[j];

        temp_counter++;
      } // end of if-statement
    } // end of inner loop
  } // end of outer loop
  row_ptr[nodes] = nz;

  free(FromNode);
  free(ToNode);
  free(num_out_links);

} // end of read_graph_file

void spmm(int nodes, double *val, int *col_idx, int *row_ptr,
          double *x_old, double *spmm_result, double d){
  /*Performing the Sparse Matrix vector multiplictation A*x^(k-1)
  spmm_result is the vector for storing the result of the multiplictation*/

  int i, j;
  for (i = 0; i < nodes; i++){
    spmm_result[i] = 0; // Initilize

    for (j = row_ptr[i]; j <= row_ptr[i+1] - 1; j++){
        spmm_result[i] += val[j] * x_old[col_idx[j]];

    } //end of inner loop

    spmm_result[i] *= d;

  }   // end of outer loop
}     // end of function

void PageRank_iterations(double *val, int *col_idx,
                         int *row_ptr, int *dangling_idx,
                         int nodes, double d, double W_sum,
                         double *x_old, double *x_new,
                         double *spmm_result, double one_over_N) {


  W_sum = dangling_PageRank(x_old, dangling_idx);
  spmm(nodes, val, col_idx, row_ptr, x_old, spmm_result, d);
  //W_sum og spmm beregnes bare én gang per iterasjon av PageRank_iterations

  for (int i = 0; i < nodes; i++){
    x_new[i] = ((1 - d + d*W_sum)*one_over_N) + spmm_result[i];
  }

  // for (int i = 0; i < nodes; i++) printf("%d  |  %lf\n", i, x_old[i]);
  // printf("\n New iteration \n");
}
