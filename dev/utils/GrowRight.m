% GrowRight.m
% contracts a given dmpo site tensor into the right block -- a
% contraction from the last site in the system through to the site after
% this one
% Oliver Thomson Brown
% 2016-03-22
%
% updateBlock = GrowRight(siteTensor, mpo, rightBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
%
% RETURN
% updateBlock       : rank-3 complex double array, contains the contraction from
%                     the last site through to this one
%
% INPUTS
% siteTensor        : rank-4 complex double array, the dmpo tensor for the
%                     specified site
% mpo               : rank-6 complex double array, mpo for the operator being
%                     evaluated in the form <rho|| O ||rho>
% rightBlock        : rank-3 complex double array, contains the contraction from
%                     the last site through to the site immediately after this
%                     one
% ROW_SIZE          : integer, the size of the first virtual dimension of
%                     siteTensor
% COL_SIZE          : integer, the size of the second virtual dimension of
%                     siteTensor
% OP_ROW            : integer, the size of the first virtual dimension of mpo
% OP_COL            : integer, the size of the second virtual dimension of mpo 

function [updateBlock] = GrowRight(siteTensor, mpo, rightBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
    updateBlock = zeros(ROW_SIZE, OP_ROW, ROW_SIZE);

    conjTensor = zeros(COL_SIZE, ROW_SIZE, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            conjTensor(:, :, bra, ket) = ctranspose(siteTensor(:, :, bra, ket));
        end
    end

    for conjCol = 1 : 1 : ROW_SIZE
        for opRow = 1 : 1 : OP_ROW
            for row = 1 : 1 : ROW_SIZE
                BWFB = 0;
                for conjBra = 1 : 1 : HILBY
                    for conjKet = 1 : 1 : HILBY
                        for conjRow = 1 : 1 : COL_SIZE
                            WFB = 0;
                            for bra = 1 : 1 : HILBY
                                for ket = 1 : 1 : HILBY
                                    for opCol = 1 : 1 : OP_COL
                                        FB = 0;
                                        for col = 1 : 1 : COL_SIZE
                                            FB = FB + rightBlock(conjRow, opCol, col) * siteTensor(row, col, bra, ket);
                                        end
                                        WFB = WFB + mpo(bra, ket, conjBra, conjKet, opRow, opCol) * FB;
                                    end
                                end
                            end
                            BWFB = BWFB + conjTensor(conjRow, conjCol, conjBra, conjKet) * WFB;
                        end
                    end
                end
                updateBlock(conjCol, opRow, row) = BWFB;
            end
        end
    end
end
