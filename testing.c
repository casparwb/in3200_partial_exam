# include <stdio.h>
# include <stdlib.h>
# include <omp.h>

void somefunction(int *z_ptr){
  *z_ptr = 3;
}

int main(){
  int z;
  int *z_ptr;
  z_ptr = &z;
  z = 0;

  somefunction(z_ptr);


  return 0;
}
