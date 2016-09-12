% EffOpTensor.m
% function which returns the effective operator tensor for some particular site
% it performs the contraction with the left block tensor, right block tensor,
% and matrix product operator tensor for the site
% Oliver Thomson Brown
% 2016-05-04
%
% [effOpTensor] = EffOpTensor(lBlock, mpo, rBlock, ROW_SIZE, COL_SIZE, HILBY)
%
% RETURN
% effOpTensor     : 8-dimensional complex double array, rank-6 effective
%                   operator for some site. Indices are effOpTensor(row, col,
%                   conjRow, conjCol, bra, ket, conjBra, conjKet)
% INPUTS
% lBlock          : 3-dimensional complex double array, rank-3 left block tensor
%                   which contains the contraction through the system from the
%                   first site up to the target site
% siteMPO         : 6-dimensional complex double array, rank-6 matrix product
%                   operator tensor for the target site
% rBlock          : 3-dimensional complex double array, rank-3 right block
%                   right block tensor which contains the contraction through
%                   the sytem from the target site to the last
% ROW_SIZE        : double, the size of the first virtual dimension of the
%                   density matrix product operator for the target site
% COL_SIZE        : double, the size of the second virtual dimension of the
%                   density matrix product operator for the target site
% HILBY           : double, the size of the local state space

function [effOpTensor] = EffOpTensor(lBlock, siteMPO, rBlock, ROW_SIZE, COL_SIZE, HILBY)
    % effOpTensor(row, col, conjRow, conjCol, conjBra, conjKet, bra, ket)
    effOpTensor = zeros(ROW_SIZE, COL_SIZE, COL_SIZE, ROW_SIZE, ...
                        HILBY, HILBY, HILBY, HILBY);

    % permute input arrays to avoid for-loops
    % original indexing:
    %   lBlock(conjCol, opRow, row)
    %   siteMPO(conjBra, conjKet, bra, ket, opRow, opCol)
    %   rBlock(conjRow, opCol, col)
    lBlock = permute(lBlock, [3, 2, 1]);
    siteMPO = permute(siteMPO, [5, 6, 1, 2, 3, 4]);
    rBlock = permute(rBlock, [2, 3, 1]);

    for conjKet = 1 : 1 : HILBY
        for conjBra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                for bra = 1 : 1 : HILBY
                    for conjRow = 1 : 1 : COL_SIZE
                        for conjCol = 1 : 1 : ROW_SIZE
                            effOpTensor(:, :, conjRow, conjCol, bra, ket, ...
                                        conjBra, conjKet) = ...
                            lBlock(:, :, conjCol) ...
                            * siteMPO(:, :, conjBra, conjKet, bra, ket) ...
                            * rBlock(:, :, conjRow);
                        end
                    end
                end
            end
        end
    end
end
