
void read_graph_file(FILE *infile, double **val_ptr, int **col_idx_ptr,
              int **row_ptr_ptr, int **dangling_idx_ptr, int z, int nodes, int edges) {

  int temp1, temp2, temp_counter;
  int i, j, k, nz = 0;  // nz: counter for no. of nonzero elements

  char buffer[100]; // Array for storing lines that we have to skip

  int *FromNode = (int*)malloc(edges*sizeof(int));      //array for storing the FromNodeId-values
  int *ToNode = (int*)malloc(edges*sizeof(int));        // array for storing the ToNodeId-values
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
    *dangling_idx_ptr = (int*)malloc(z*sizeof(int));
    temp_counter = 0;
    for (i = 0; i < nodes; i++){
      if (num_out_links[i] == 0){
        (*dangling_idx_ptr)[temp_counter] = i;
        temp_counter++;
      } // end of inner if-statement
    }   // end of inner loop
  }     // end of outer if-statement

  /* Building the hyperlink matrix in CRS-format */
  *val_ptr     = (double*) malloc(nz*sizeof(double));
  *col_idx_ptr = (int*)    malloc(nz*sizeof(int));
  *row_ptr_ptr = (int*)    malloc((nodes+1)*sizeof(int));


  temp_counter = 0;
  for (i = 0; i < nodes; i++){
    (*row_ptr_ptr)[i] = temp_counter;

    for (j = 0; j < edges; j++){

      if ((ToNode[j] == i) && (ToNode[j] != FromNode[j])){
        (*val_ptr)[temp_counter] = 1.0/num_out_links[FromNode[j]];
        (*col_idx_ptr)[temp_counter] = FromNode[j];

        temp_counter++;
      } // end of if-statement
    } // end of inner loop
  } // end of outer loop
  (*row_ptr_ptr)[nodes] = nz;

  free(FromNode);
  free(ToNode);
  free(num_out_links);

} // end of read_graph_file


void PageRank_iterations(double *val, int *col_idx,
                         int *row_ptr, int *dangling_idx,
                         int nodes, double d, double W_sum,
                         double *x_old, double *x_new,
                         double *spmm_result, double one_over_N, int z) {




    W_sum = 0;
    if (z = 0) dangling_idx[0] = 0; // to stop segmentation fault
    else{
      #pragma omp parallel for reduction(+:W_sum)
      for (int i = 0; i < z; i++) {
        W_sum += x_old[dangling_idx[i]];
      }
    }

    for (int i = 0; i < nodes; i++){
      spmm_result[i] = 0; // Initilize

      for (int j = row_ptr[i]; j < row_ptr[i+1]; j++){
          spmm_result[i] += val[j] * x_old[col_idx[j]];
      } //end of inner loop

      spmm_result[i] *= d;
    }

    for (int i = 0; i < nodes; i++){
      x_new[i] = ((1 - d + d*W_sum)*one_over_N) + spmm_result[i];
    }
} // end of PageRank-function

void top_n_pages(double *x, double *top_pages, int *top_pages_idx, int n, int nodes){
  /*Finding the top n PageRank scores and their corresponding website indices*/
  double temp_max, max_idx, max = 1;

  for (int j = 0; j < n; j++){
    temp_max = 0;

    for (int i = 0; i < nodes; i++){
      if ((x[i] > temp_max) && (x[i] < max)){
        temp_max = x[i];
        max_idx = i;
      } // end of if-statement
    } // end of inner for loop

    max = temp_max;
    top_pages_idx[j] = max_idx;
    top_pages[j] = max;

  } // end of outer for loop

} // end of function
