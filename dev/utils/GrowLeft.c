/*
  GrowLeft.c
  contracts a given dmpo site tensor into the left block, which is the
  contraction of the tensor network from the first site through to the site
  prior to this one.
  Oliver Thomson Brown
  2017-03-03
*/

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <complex.h>

typedef struct array_Z3D
{
  mwSize dims[3];
  complex double ***arr;
} az3_t;

typedef struct array_Z4D
{
  mwSize dims[4];
  complex double ****arr;
} az4_t;

typedef struct array_Z6D
{
  mwSize dims[6];
  complex double ******arr;
} az6_t;

complex double ***init_Z3arr(mwSize dims[3])
{
  mwIndex dim0, dim1, dim2;

  complex double ***arr = mxMalloc(sizeof(complex double **) * dims[0]);
  for(dim0 = 0; dim0 < dims[0]; ++dim0) {
    arr[dim0] = mxMalloc(sizeof(complex double *) * dims[1]);
    for(dim1 = 0; dim1 < dims[1]; ++dim1) {
      arr[dim0][dim1] = mxMalloc(sizeof(complex double) * dims[2]);
    }
  }

  return arr;
}

complex double ****init_Z4arr(mwSize dims[4])
{
  mwIndex dim0, dim1, dim2;

  complex double ****arr = mxMalloc(sizeof(complex double ***) * dims[0]);
  for(dim0 = 0; dim0 < dims[0]; ++dim0) {
    arr[dim0] = mxMalloc(sizeof(complex double **) * dims[1]);
    for(dim1 = 0; dim1 < dims[1]; ++dim1) {
      arr[dim0][dim1] = mxMalloc(sizeof(complex double *) * dims[2]);
      for(dim2 = 0; dim2 < dims[2]; ++dim2) {
        arr[dim0][dim1][dim2] = mxMalloc(sizeof(complex double) * dims[3]);
      }
    }
  }
  return arr;
}

complex double ******init_Z6arr(mwSize dims[6])
{
  mwIndex dim0, dim1, dim2, dim3, dim4;

  complex double ******arr = mxMalloc(sizeof(complex double *****) * dims[0]);
  for(dim0 = 0; dim0 < dims[0]; ++dim0) {
    arr[dim0] = mxMalloc(sizeof(complex double ****) * dims[1]);
    for(dim1 = 0; dim1 < dims[1]; ++dim1) {
      arr[dim0][dim1] = mxMalloc(sizeof(complex double ***) * dims[2]);
      for(dim2 = 0; dim2 < dims[2]; ++dim2) {
        arr[dim0][dim1][dim2] = mxMalloc(sizeof(complex double **) * dims[3]);
        for(dim3 = 0; dim3 < dims[3]; ++dim3) {
          arr[dim0][dim1][dim2][dim3]
                                = mxMalloc(sizeof(complex double *) * dims[4]);
          for(dim4 = 0; dim4 < dims[4]; ++dim4) {
            arr[dim0][dim1][dim2][dim3][dim4]
                                = mxMalloc(sizeof(complex double) * dims[5]);
          }
        }
      }
    }
  }

  return arr;
}

void free_Z3arr(az3_t *X)
{
  mwIndex dim, dim0, dim1;

  mwSize dims[4];
  for(dim = 0; dim < 4; ++dim) {
    dims[dim] = (*X).dims[dim];
  }

  for(dim0 = (dims[0] - 1); dim0 > -1; --dim0) {
    for(dim1 = (dims[1] - 1); dim1 > -1; --dim1) {
      mxFree((*X).arr[dim0][dim1]);
    }
    mxFree((*X).arr[dim0]);
  }
  mxFree((*X).arr);

  return;
}

void free_Z4arr(az4_t *X)
{
  mwIndex dim, dim0, dim1, dim2;

  mwSize dims[4];
  for(dim = 0; dim < 4; ++dim) {
    dims[dim] = (*X).dims[dim];
  }

  for(dim0 = (dims[0] - 1); dim0 > -1; --dim0) {
    for(dim1 = (dims[1] - 1); dim1 > -1; --dim1) {
      for(dim2 = (dims[2] - 1); dim2 > -1; --dim2) {
        mxFree((*X).arr[dim0][dim1][dim2]);
      }
      mxFree((*X).arr[dim0][dim1]);
    }
    mxFree((*X).arr[dim0]);
  }
  mxFree((*X).arr);

  return;
}

void free_Z6arr(az6_t *X)
{
  mwIndex dim, dim0, dim1, dim2, dim3, dim4;

  mwSize dims[6];
  for(dim = 0; dim < 6; ++dim) {
    dims[dim] = (*X).dims[dim];
  }

  for(dim0 = (dims[0] - 1); dim0 > -1; --dim0) {
    for(dim1 = (dims[1] - 1); dim1 > -1; --dim1) {
      for(dim2 = (dims[2] - 1); dim2 > -1; --dim2) {
        for(dim3 = (dims[3] - 1); dim3 > -1; --dim3) {
          for(dim4 = (dims[4] - 1); dim4 > -1; --dim4) {
            mxFree((*X).arr[dim0][dim1][dim2][dim3][dim4]);
          }
          mxFree((*X).arr[dim0][dim1][dim2][dim3]);
        }
        mxFree((*X).arr[dim0][dim1][dim2]);
      }
      mxFree((*X).arr[dim0][dim1]);
    }
    mxFree((*X).arr[dim0]);
  }
  mxFree((*X).arr);

  return;
}

az3_t *GrowLeft(const az4_t *siteTensor, az6_t *mpo, const az3_t *leftBlock,
                const mwSize ROW_SIZE, const mwSize COL_SIZE,
                const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL)
{
  mwIndex row, col, conjRow, conjCol, bra, ket;
  mwIndex conjBra, conjKet, opRow, opCol, dim;
  complex double FA, WFA, AWFA;

  /* create conjugate site tensor */
  az4_t *conjTensor = mxMalloc(sizeof(az4_t));
  for(dim = 0; dim < 4; ++dim) {
    (*conjTensor).dims[dim] = (*siteTensor).dims[dim];
  }
  (*conjTensor).arr = init_Z4arr((*conjTensor).dims);
  for(row = 0; row < ROW_SIZE; ++row) {
    for(col = 0; col < COL_SIZE; ++col) {
      for(bra = 0; bra < HILBY; ++bra) {
        for(ket = 0; ket < HILBY; ++ket) {
          (*conjTensor).arr[row][col][bra][ket]
            = conj((*siteTensor).arr[row][col][bra][ket]);
        }
      }
    }
  }

  /* allocate return */
  az3_t *updateBlock = mxMalloc(sizeof(az3_t));
  (*updateBlock).dims[0] = COL_SIZE;
  (*updateBlock).dims[1] = OP_COL;
  (*updateBlock).dims[2] = COL_SIZE;
  (*updateBlock).arr = init_Z3arr((*updateBlock).dims);

  free_Z4arr(conjTensor);

  return updateBlock;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return;
}
