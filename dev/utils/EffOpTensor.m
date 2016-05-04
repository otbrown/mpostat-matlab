% EffOpTensor.m
% function which returns the effective operator tensor for some particular site
% it performs the contraction with the left block tensor, right block tensor,
% and matrix product operator tensor for the site
% Oliver Thomson Brown
% 2016-05-04
%
% [effOpTensor] = EffOpTensor(lBlock, mpo, rBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
%
% RETURN
% effOpTensor     : 8-dimensional complex double array, rank-6 effective
%                   operator for some site. Indices are effOpTensor(row, col,
%                   conjRow, conjCol, bra, ket, conjBra, conjKet)
% INPUTS
%

function [effOpTensor] = EffOpTensor(lBlock, mpo, rBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
    effOpTensor = zeros(ROW_SIZE, COL_SIZE, COL_SIZE, ROW_SIZE, ...
                        HILBY, HILBY, HILBY, HILBY);

    % naieve for looping
    for conjKet = 1 : 1 : HILBY
        for conjBra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                for bra = 1 : 1 : HILBY
                    for conjCol = 1 : 1 : ROW_SIZE
                        for conjRow = 1 : 1 : COL_SIZE
                            for col = 1 : 1 : COL_SIZE
                                for row = 1 : 1 : ROW_SIZE
                                    LOR = 0;
                                    for opRow = 1 : 1 : OP_ROW
                                        for opCol = 1 : 1 : OP_COL
                                            LOR = LOR ...
                                            + lBlock(row, opRow, conjCol) ...
                                            * mpo(conjBra, conjKet, ...
                                            bra, ket, opRow, opCol) ...
                                            * rBlock(col, opCol, conjRow);
                                        end
                                    end
                                    effOpTensor(row, col, conjRow, conjCol, ...
                                    bra, ket, conjBra, conjKet) = LOR;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
