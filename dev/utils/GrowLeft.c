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
  complex double *arr;
} az3_t;

typedef struct array_Z4D
{
  mwSize dims[4];
  complex double *arr;
} az4_t;

typedef struct array_Z6D
{
  mwSize dims[6];
  complex double *arr;
} az6_t;

az3_t *init_Z3arr(const mwSize dims[3])
{
  const mwSize NUMEL = dims[0] * dims[1] * dims[2];
  az3_t *A = mxMalloc(sizeof(az3_t));
  mwIndex dim;

  for(dim = 0; dim < 3; ++dim) {
    (*A).dims[dim] = dims[dim];
  }

  (*A).arr = mxMalloc(sizeof(complex double) * NUMEL);

  return A;
}

az4_t *init_Z4arr(const mwSize dims[4])
{
  const mwSize NUMEL = dims[0] * dims[1] * dims[2] * dims[3];
  az4_t *A = mxMalloc(sizeof(az4_t));
  mwIndex dim;

  for(dim = 0; dim < 4; ++dim) {
    (*A).dims[dim] = dims[dim];
  }

  (*A).arr = mxMalloc(sizeof(complex double) * NUMEL);

  return A;
}

az6_t *init_Z6arr(const mwSize dims[6])
{
  const mwSize NUMEL = dims[0] * dims[1] * dims[2]
                        * dims[3] * dims[4] * dims[5];
  az6_t *A = mxMalloc(sizeof(az6_t));
  mwIndex dim;

  for(dim = 0; dim < 6; ++dim) {
    (*A).dims[dim] = dims[dim];
  }

  (*A).arr = mxMalloc(sizeof(complex double) * NUMEL);

  return A;
}

void free_Z3arr(az3_t *X)
{
  mxFree((*X).arr);
  mxFree(X);

  return;
}

void free_Z4arr(az4_t *X)
{
  mxFree((*X).arr);
  mxFree(X);

  return;
}

void free_Z6arr(az6_t *X)
{
  mxFree((*X).arr);
  mxFree(X);

  return;
}

az3_t *MxToZ3Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const mwSize NDIMS = mxGetNumberOfDimensions(M);
  const bool isComplex = mxIsComplex(M);
  double *MXDATAr, *MXDATAi;
  az3_t *A;
  mwIndex lindex, dim;
  mwSize dims[3];

  if(isComplex) {
    MXDATAr = mxGetPr(M);
    MXDATAi = mxGetPi(M);
  } else {
    MXDATAr = mxGetPr(M);
  }

  /* pad any missing dimensions with 1 */
  if(NDIMS < 3) {
    for(dim = 0; dim < NDIMS; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
    for(dim = NDIMS; dim < 3; ++dim) {
      dims[dim] = 1;
    }
  }
  else {
    for(dim = 0; dim < 3; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
  }
  const mwSize NUMEL = dims[0] * dims[1] * dims[2];

  /* initialise return struct */
  A = init_Z3arr(dims);

  /* fill from mxArray, converting to complex double */
  if(isComplex) {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * MXDATAi[lindex];
    }
  } else {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * 0.0;
    }
  }

  return A;
}

az4_t *MxToZ4Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const mwSize NDIMS = mxGetNumberOfDimensions(M);
  const bool isComplex = mxIsComplex(M);
  double *MXDATAr, *MXDATAi;
  az4_t *A;
  mwIndex lindex, dim;
  mwSize dims[4];

  if(isComplex) {
    MXDATAr = mxGetPr(M);
    MXDATAi = mxGetPi(M);
  } else {
    MXDATAr = mxGetPr(M);
  }

  /* pad any missing dimensions with 1 */
  if(NDIMS < 4) {
    for(dim = 0; dim < NDIMS; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
    for(dim = NDIMS; dim < 4; ++dim) {
      dims[dim] = 1;
    }
  } else {
    for(dim = 0; dim < 4; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
  }
  const mwSize NUMEL = dims[0] * dims[1] * dims[2] * dims[3];

  /* initialise return struct */
  A = init_Z4arr(dims);
  /* fill from mxArray, converting to complex double */
  if(isComplex) {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * MXDATAi[lindex];
    }
  } else {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * 0.0;
    }
  }

  return A;
}

az6_t *MxToZ6Array(const mxArray *M)
{
  const mwSize *MXDIMS = mxGetDimensions(M);
  const mwSize NDIMS = mxGetNumberOfDimensions(M);
  const bool isComplex = mxIsComplex(M);
  double *MXDATAr, *MXDATAi;
  az6_t *A;
  mwIndex lindex, dim;
  mwSize dims[6];

  if(isComplex) {
    MXDATAr = mxGetPr(M);
    MXDATAi = mxGetPi(M);
  } else {
    MXDATAr = mxGetPr(M);
  }

  /* pad any missing dimensions with 1 */
  if(NDIMS < 6) {
    for(dim = 0; dim < NDIMS; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
    for(dim = NDIMS; dim < 6; ++dim) {
      dims[dim] = 1;
    }
  } else {
    for(dim = 0; dim < 6; ++dim) {
      dims[dim] = MXDIMS[dim];
    }
  }
  const mwSize NUMEL = dims[0] * dims[1] * dims[2]
                        * dims[3] * dims[4] * dims[5];

  /* initialise return struct */
  A = init_Z6arr(dims);

  /* fill from mxArray, converting to complex double */
  if(isComplex) {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * MXDATAi[lindex];
    }
  } else {
    for(lindex = 0; lindex < NUMEL; ++lindex) {
      (*A).arr[lindex] = (complex double) MXDATAr[lindex] + I * 0.0;
    }
  }

  return A;
}

mxArray *Z3ArrayToMx(az3_t *A)
{
  mxArray *M = mxCreateUninitNumericArray(3, (*A).dims, mxDOUBLE_CLASS,
                                          mxCOMPLEX);
  const mwSize NUMEL = (*A).dims[0] * (*A).dims[1] * (*A).dims[2];
  mwIndex lindex;

  double *mxVecR = mxMalloc(sizeof(double) * NUMEL);
  double *mxVecI = mxMalloc(sizeof(double) * NUMEL);

  for(lindex = 0; lindex < NUMEL; ++lindex) {
    mxVecR[lindex] = creal((*A).arr[lindex]);
    mxVecI[lindex] = cimag((*A).arr[lindex]);
  }

  mxSetPr(M, mxVecR);
  mxSetPi(M, mxVecI);

  return M;
}

mwIndex a3dex(const mwIndex dim0, const mwIndex dim1, const mwIndex dim2,
              const mwSize dims[3])
{
  mwIndex lindex = dim2 * dims[1] * dims[0] + dim1 * dims[0] + dim0;

  return lindex;
}

mwIndex a4dex(const mwIndex dim0, const mwIndex dim1, const mwIndex dim2,
              const mwIndex dim3, const mwSize dims[4])
{
  mwIndex lindex = dim3 * dims[2] * dims[1] * dims[0] + dim2 * dims[1] * dims[0]
                   + dim1 * dims[0] + dim0;

  return lindex;
}

mwIndex a6dex(const mwIndex dim0, const mwIndex dim1, const mwIndex dim2,
              const mwIndex dim3, const mwIndex dim4, const mwIndex dim5,
              const mwSize dims[6])
{
  mwIndex lindex = dim5 * dims[4] * dims[3] * dims[2] * dims[1] * dims[0]
                   + dim4 * dims[3] * dims[2] * dims[1] * dims[0]
                   + dim3 * dims[2] * dims[1] * dims[0]
                   + dim2 * dims[1] * dims[0] + dim1 * dims[0] + dim0;

  return lindex;
}

az3_t *GrowLeft(az4_t *siteTensor, az6_t *mpo, az3_t *leftBlock,
                const mwSize ROW_SIZE, const mwSize COL_SIZE,
                const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL)
{
  mwIndex row, col, conjRow, conjCol, bra, ket;
  mwIndex conjBra, conjKet, opRow, opCol;
  mwIndex sitedex, conjdex, mpodex, leftdex, updex;
  mwSize numel, updims[3];
  complex double FA, WFA, AWFA;
  az3_t *updateBlock;

  /* create conjugate site tensor */
  az4_t *conjTensor = init_Z4arr((*siteTensor).dims);
  numel = (*siteTensor).dims[0] * (*siteTensor).dims[1]
          * (*siteTensor).dims[2] * (*siteTensor).dims[3];
  for(conjdex = 0; conjdex < numel; ++conjdex) {
    (*conjTensor).arr[conjdex] = conj((*siteTensor).arr[conjdex]);
  }

  /* allocate return */
  updims[0] = COL_SIZE;
  updims[1] = OP_COL;
  updims[2] = COL_SIZE;
  updateBlock = init_Z3arr(updims);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFA = 0.0 + I * 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, (*conjTensor).dims);

              WFA = 0.0 + I * 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, (*mpo).dims);

                    FA = 0.0 + I * 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, (*leftBlock).dims);
                      sitedex = a4dex(row, col, bra, ket, (*siteTensor).dims);

                      FA += (*leftBlock).arr[leftdex] * (*siteTensor).arr[sitedex];
                    }
                    WFA += (*mpo).arr[mpodex] * FA;
                  }
                }
              }
              AWFA += (*conjTensor).arr[conjdex] * WFA;
            }
          }
        }
        (*updateBlock).arr[updex] = AWFA;
        ++updex;
      }
    }
  }

  free_Z4arr(conjTensor);

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

  return;
}
