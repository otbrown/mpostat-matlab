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

complex double ***init_Z3arr(const mwSize dims[3])
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

complex double ****init_Z4arr(const mwSize dims[4])
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

complex double ******init_Z6arr(const mwSize dims[6])
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

az3_t *MxToZ3Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const double *MXDATAr = mxGetPr(M);
  const double *MXDATAi = mxGetPi(M);
  mwIndex dim, dim0, dim1, dim2, lindex;

  /* initialise return struct */
  az3_t *A = mxMalloc(sizeof(az3_t));
  for(dim = 0; dim < 3; ++dim) {
    (*A).dims[dim] = MXDIMS[dim];
  }
  (*A).arr = init_Z3arr(MXDIMS);

  lindex = 0;
  for(dim2 = 0; dim2 < (*A).dims[2]; ++dim2) {
    for(dim1 = 0; dim1 < (*A).dims[1]; ++dim1) {
      for(dim0 = 0; dim0 < (*A).dims[0]; ++dim0) {
        (*A).arr[dim0][dim1][dim2] = MXDATAr[lindex] + I * MXDATAi[lindex];
        ++lindex;
      }
    }
  }

  return A;
}

az4_t *MxToZ4Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const double *MXDATAr = mxGetPr(M);
  const double *MXDATAi = mxGetPi(M);
  mwIndex dim, dim0, dim1, dim2, dim3, lindex;

  /* initialise return struct */
  az4_t *A = mxMalloc(sizeof(az4_t));
  for(dim = 0; dim < 4; ++dim) {
    (*A).dims[dim] = MXDIMS[dim];
  }
  (*A).arr = init_Z4arr(MXDIMS);

  lindex = 0;
  for(dim3 = 0; dim3 < (*A).dims[3]; ++dim3) {
    for(dim2 = 0; dim2 < (*A).dims[2]; ++dim2) {
      for(dim1 = 0; dim1 < (*A).dims[1]; ++dim1) {
        for(dim0 = 0; dim0 < (*A).dims[0]; ++dim0) {
          (*A).arr[dim0][dim1][dim2][dim3]
                                        = MXDATAr[lindex] + I * MXDATAi[lindex];
          ++lindex;
        }
      }
    }
  }

  return A;
}

az6_t *MxToZ6Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const double *MXDATAr = mxGetPr(M);
  const double *MXDATAi = mxGetPi(M);
  mwIndex dim, dim0, dim1, dim2, dim3, dim4, dim5, lindex;

  /* initialise return struct */
  az6_t *A = mxMalloc(sizeof(az6_t));
  for(dim = 0; dim < 6; ++dim) {
    (*A).dims[dim] = MXDIMS[dim];
  }
  (*A).arr = init_Z6arr(MXDIMS);

  lindex = 0;
  for(dim5 = 0; dim5 < (*A).dims[5]; ++dim5) {
    for(dim4 = 0; dim4 < (*A).dims[4]; ++dim4) {
      for(dim3 = 0; dim3 < (*A).dims[3]; ++dim3) {
        for(dim2 = 0; dim2 < (*A).dims[2]; ++dim2) {
          for(dim1 = 0; dim1 < (*A).dims[1]; ++dim1) {
            for(dim0 = 0; dim0 < (*A).dims[0]; ++dim0) {
              (*A).arr[dim0][dim1][dim2][dim3][dim4][dim5]
                                        = MXDATAr[lindex] + I * MXDATAi[lindex];
              ++lindex;
            }
          }
        }
      }
    }
  }

  return A;
}

mxArray *Z3ArrayToMx(az3_t *A)
{
  mxArray *M = mxCreateUninitNumericArray(3, (*A).dims, mxDOUBLE_CLASS,
                                          mxCOMPLEX);
  const mwSize NUMEL = (*A).dims[0] * (*A).dims[1] * (*A).dims[2];
  mwIndex dim0, dim1, dim2, lindex;

  double *mxVecR = mxMalloc(sizeof(double) * NUMEL);
  double *mxVecI = mxMalloc(sizeof(double) * NUMEL);

  lindex = 0;
  for(dim2 = 0; dim2 < (*A).dims[2]; ++dim2) {
    for(dim1 = 0; dim1 < (*A).dims[1]; ++dim1) {
      for(dim0 = 0; dim0 < (*A).dims[0]; ++dim0) {
        mxVecR[lindex] = creal((*A).arr[dim0][dim1][dim2]);
        mxVecI[lindex] = cimag((*A).arr[dim0][dim1][dim2]);
        ++lindex;
      }
    }
  }

  mxSetPr(M, mxVecR);
  mxSetPi(M, mxVecI);

  return M;
}

az3_t *GrowLeft(az4_t *siteTensor, az6_t *mpo, az3_t *leftBlock,
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

  /* loop the loop! */
  for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(col = 0; col < COL_SIZE; ++col) {
        AWFA = 0.0 + I * 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              WFA = 0.0 + I * 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    FA = 0.0 + I * 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      FA += (*leftBlock).arr[conjCol][opRow][row] * (*siteTensor).arr[row][col][bra][ket];
                    }
                    WFA += (*mpo).arr[conjBra][conjKet][bra][ket][opRow][opCol] * FA;
                  }
                }
              }
              AWFA += (*conjTensor).arr[conjCol][conjRow][conjBra][conjKet] * WFA;
            }
          }
        }
        (*updateBlock).arr[conjRow][opCol][col] = AWFA;
      }
    }
  }

  free_Z4arr(conjTensor);
  mxFree(conjTensor);

  return updateBlock;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  az3_t *updateBlock;
  az4_t *siteTensor;
  az6_t *mpo;
  az3_t *leftBlock;

  siteTensor = MxToZ4Array(prhs[0]);
  mpo = MxToZ6Array(prhs[1]);
  leftBlock = MxToZ3Array(prhs[2]);

  const mwSize ROW_SIZE = (*siteTensor).dims[0];
  const mwSize COL_SIZE = (*siteTensor).dims[1];
  const mwSize HILBY = (*siteTensor).dims[2];
  const mwSize OP_ROW = (*mpo).dims[4];
  const mwSize OP_COL = (*mpo).dims[5];

  updateBlock = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE,
                          HILBY, OP_ROW, OP_COL);

  plhs[0] = Z3ArrayToMx(updateBlock);

  free_Z3arr(updateBlock);
  free_Z4arr(siteTensor);
  free_Z6arr(mpo);
  free_Z3arr(leftBlock);
  mxFree(updateBlock);
  mxFree(siteTensor);
  mxFree(mpo);
  mxFree(leftBlock);

  return;
}
