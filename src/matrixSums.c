#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

SEXP _colSumByGroup(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int *x = INTEGER(R_x);

  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group)) {
    error("The grouping argument must be a factor");
  }

  int *group = INTEGER(R_group);
  int nl = nlevels(R_group);
  // If the sizes of the grouping variable and matrix do not match, throw an error
  if (LENGTH(R_group) != nc) {
    error("The length of the grouping argument must match the number of columns in the matrix");
  }
  for (i = 0; i < nc; i++) {
    if(group[i] == NA_INTEGER) {
      error("Labels in group and pgroup must not be NA.");
    }
  }  
  
  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(INTSXP, nr, nl));
  // Set a pointer to the return matrix and initialize the memory
  int *ans = INTEGER(R_ans);
  Memzero(ans, nr * nl);

  // Sum the totals for each element of the 'group' variable
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) {
      ans[(group[j] - 1) * nr + i] += x[j * nr + i];
    }
  }

  UNPROTECT(1);
  return(R_ans);
}

SEXP _colSumByGroup_numeric(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  double *x = REAL(R_x);

  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group)) {
    error("The grouping argument must be a factor");
  }

  int *group = INTEGER(R_group);
  int nl = nlevels(R_group);
  // If the sizes of the grouping variable and matrix do not match, throw an error
  if (LENGTH(R_group) != nc) {
    error("The length of the grouping argument must match the number of columns in the matrix");
  }
  
  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(REALSXP, nr, nl));
  // Set a pointer to the return matrix and initialize the memory
  double *ans = REAL(R_ans);
  Memzero(ans, nr * nl);

  // Sum the totals for each element of the 'group' variable
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) {
      ans[(group[j] - 1) * nr + i] += x[j * nr + i];
    }
  }

  UNPROTECT(1);
  return(R_ans);
}
