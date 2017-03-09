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

double zmultr(const double Ar, const double Ai, const double Br, const double Bi)
{
  double Zr = (Ar * Br) - (Ai * Bi);

  return Zr;
}

double zmulti(const double Ar, const double Ai, const double Br, const double Bi)
{
  double Zi = (Ar * Bi) + (Ai * Br);

  return Zi;
}

mxArray * GrowLeft_RRR(const double* siteReal, const double* mpoReal, const double* leftReal, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxREAL);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, WFAr, AWFAr;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                    }
                    WFAr += mpoReal[mpodex] * FAr;
                  }
                }
              }
              AWFAr += siteReal[conjdex] * WFAr;
            }
          }
        }
        updReal[updex] = AWFAr;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);

  return updateBlock;
}

mxArray *GrowLeft_RRZ(const double* siteReal, const double* mpoReal, const double* leftReal, const double* leftImag, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                      FAi += leftImag[leftdex] * siteReal[sitedex];
                    }
                    WFAr += mpoReal[mpodex] * FAr;
                    WFAi += mpoReal[mpodex] * FAi;
                  }
                }
              }
              AWFAr += siteReal[conjdex] * WFAr;
              AWFAi += siteReal[conjdex] * WFAi;
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_RZR(const double* siteReal, const double* mpoReal, const double* mpoImag, const double* leftReal, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                    }
                    WFAr += mpoReal[mpodex] * FAr;
                    WFAi += mpoImag[mpodex] * FAr;
                  }
                }
              }
              AWFAr += siteReal[conjdex] * WFAr;
              AWFAi += siteReal[conjdex] * WFAi;
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_RZZ(const double* siteReal, const double* mpoReal, const double* mpoImag, const double* leftReal, const double* leftImag, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                      FAi += leftImag[leftdex] * siteReal[sitedex];
                    }
                    WFAr += zmultr(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                    WFAi += zmulti(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                  }
                }
              }
              AWFAr += siteReal[conjdex] * WFAr;
              AWFAi += siteReal[conjdex] * WFAi;
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_ZRR(const double* siteReal, const double* siteImag, const double* mpoReal, const double* leftReal, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                      FAi += leftReal[leftdex] * siteImag[sitedex];
                    }
                    WFAr += mpoReal[mpodex] * FAr;
                    WFAi += mpoReal[mpodex] * FAi;
                  }
                }
              }
              AWFAr += zmultr(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
              AWFAi += zmulti(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_ZRZ(const double* siteReal, const double* siteImag, const double* mpoReal, const double* leftReal, const double* leftImag, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += zmultr(leftReal[leftdex], leftImag[leftdex], siteReal[sitedex], siteImag[sitedex]);
                      FAi += zmulti(leftReal[leftdex], leftImag[leftdex], siteReal[sitedex], siteImag[sitedex]);
                    }
                    WFAr += mpoReal[mpodex] * FAr;
                    WFAi += mpoReal[mpodex] * FAi;
                  }
                }
              }
              AWFAr += zmultr(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
              AWFAi += zmulti(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_ZZR(const double* siteReal, const double* siteImag, const double* mpoReal, const double* mpoImag, const double* leftReal, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += leftReal[leftdex] * siteReal[sitedex];
                      FAi += leftReal[leftdex] * siteImag[sitedex];
                    }
                    WFAr += zmultr(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                    WFAi += zmulti(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                  }
                }
              }
              AWFAr += zmultr(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
              AWFAi += zmulti(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

mxArray *GrowLeft_ZZZ(const double* siteReal, const double* siteImag, const double* mpoReal, const double* mpoImag, const double* leftReal, const double* leftImag, const mwSize ROW_SIZE, const mwSize COL_SIZE, const mwSize HILBY, const mwSize OP_ROW, const mwSize OP_COL, const mwSize* siteDims, const mwSize* mpoDims, const mwSize* lBlockDims)
{
  mwIndex row, col, conjRow, conjCol, bra, ket, conjBra, conjKet, opRow, opCol;
  mwIndex updex, sitedex, conjdex, mpodex, leftdex;
  mwSize updims[3] = {COL_SIZE, OP_COL, COL_SIZE};
  mxArray *updateBlock = mxCreateUninitNumericArray(3, updims, mxDOUBLE_CLASS, mxCOMPLEX);
  const mwSize NUMEL = COL_SIZE * OP_COL * COL_SIZE;
  double FAr, FAi, WFAr, WFAi, AWFAr, AWFAi;
  double *updReal = mxMalloc(sizeof(double) * NUMEL);
  double *updImag = mxMalloc(sizeof(double) * NUMEL);

  /* loop the loop! */
  updex = 0;
  for(col = 0; col < COL_SIZE; ++col) {
    for(opCol = 0; opCol < OP_COL; ++opCol) {
      for(conjRow = 0; conjRow < COL_SIZE; ++conjRow) {
        AWFAr = 0.0;
        AWFAi = 0.0;
        for(conjCol = 0; conjCol < ROW_SIZE; ++conjCol) {
          for(conjBra = 0; conjBra < HILBY; ++conjBra) {
            for(conjKet = 0; conjKet < HILBY; ++conjKet) {
              conjdex = a4dex(conjCol, conjRow, conjBra, conjKet, siteDims);
              WFAr = 0.0;
              WFAi = 0.0;
              for(bra = 0; bra < HILBY; ++bra) {
                for(ket = 0; ket < HILBY; ++ket) {
                  for(opRow = 0; opRow < OP_ROW; ++opRow) {
                    mpodex = a6dex(conjBra, conjKet, bra, ket, opRow, opCol, mpoDims);
                    FAr = 0.0;
                    FAi = 0.0;
                    for(row = 0; row < ROW_SIZE; ++row) {
                      leftdex = a3dex(conjCol, opRow, row, lBlockDims);
                      sitedex = a4dex(row, col, bra, ket, siteDims);
                      FAr += zmultr(leftReal[leftdex], leftImag[leftdex], siteReal[sitedex], siteImag[sitedex]);
                      FAi += zmulti(leftReal[leftdex], leftImag[leftdex], siteReal[sitedex], siteImag[sitedex]);
                    }
                    WFAr += zmultr(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                    WFAi += zmulti(mpoReal[mpodex], mpoImag[mpodex], FAr, FAi);
                  }
                }
              }
              AWFAr += zmultr(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
              AWFAi += zmulti(siteReal[conjdex], -siteImag[conjdex], WFAr, WFAi);
            }
          }
        }
        updReal[updex] = AWFAr;
        updImag[updex] = AWFAi;
        ++updex;
      }
    }
  }

  mxSetPr(updateBlock, updReal);
  mxSetPi(updateBlock, updImag);

  return updateBlock;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *siteTensor = prhs[0];
  const mxArray *mpo = prhs[1];
  const mxArray *leftBlock = prhs[2];
  const mwSize ROW_SIZE = (mwSize) mxGetScalar(prhs[3]);
  const mwSize COL_SIZE = (mwSize) mxGetScalar(prhs[4]);
  const mwSize HILBY = (mwSize) mxGetScalar(prhs[5]);
  const mwSize OP_ROW = (mwSize) mxGetScalar(prhs[6]);
  const mwSize OP_COL = (mwSize) mxGetScalar(prhs[7]);
  double *siteReal, *siteImag, *mpoReal, *mpoImag, *leftReal, *leftImag;
  mwSize numel;

  const mwSize siteDims[4] = {ROW_SIZE, COL_SIZE, HILBY, HILBY};
  const mwSize mpoDims[6] = {HILBY, HILBY, HILBY, HILBY, OP_ROW, OP_COL};
  const mwSize leftDims[3] = {ROW_SIZE, OP_ROW, ROW_SIZE};

  const bool siteCOMPLEX = mxIsComplex(siteTensor);
  const bool mpoCOMPLEX = mxIsComplex(mpo);
  const bool leftCOMPLEX = mxIsComplex(leftBlock);

  /* gather real arrays */
  siteReal = mxGetPr(siteTensor);
  mpoReal = mxGetPr(mpo);
  leftReal = mxGetPr(leftBlock);

  /* gather complex arrays, or create temporary zero arrays.. */
  if(siteCOMPLEX) {
    siteImag = mxGetPi(siteTensor);
  }

  if(mpoCOMPLEX) {
    mpoImag = mxGetPi(mpo);
  }

  if(leftCOMPLEX) {
    leftImag = mxGetPi(leftBlock);
  }

  /* run GrowLeft! */
  if(siteCOMPLEX) {
    if(leftCOMPLEX) {
      if(mpoCOMPLEX) {
        plhs[0] = GrowLeft_ZZZ(siteReal, siteImag, mpoReal, mpoImag, leftReal, leftImag, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      } else {
        plhs[0] = GrowLeft_ZRZ(siteReal, siteImag, mpoReal, leftReal, leftImag, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      }
    } else {
      if(mpoCOMPLEX) {
        plhs[0] = GrowLeft_ZZR(siteReal, siteImag, mpoReal, mpoImag, leftReal, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      } else {
        plhs[0] = GrowLeft_ZRR(siteReal, siteImag, mpoReal, leftReal, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      }
    }
  }
  else {
    if(leftCOMPLEX) {
      if(mpoCOMPLEX) {
        plhs[0] = GrowLeft_RZZ(siteReal, mpoReal, mpoImag, leftReal, leftImag, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      } else {
        plhs[0] = GrowLeft_RRZ(siteReal, mpoReal, leftReal, leftImag, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      }
    } else {
      if(mpoCOMPLEX) {
        plhs[0] = GrowLeft_RZR(siteReal, mpoReal, mpoImag, leftReal, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      } else {
        plhs[0] = GrowLeft_RRR(siteReal, mpoReal, leftReal, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL, siteDims, mpoDims, leftDims);
      }
    }
  }

  return;
}
