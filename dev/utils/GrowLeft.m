% GrowLeft.m
% contracts a given dmpo site tensor into the left block -- a contraction
% from the first site through to the site prior to this one.
% Oliver Thomson Brown
% 2016-03-22
%
% updateBlock = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE, HILBY, OP_COL)
%
% RETURN
% updateBlock   : complex double array, rank-3, contraction from the first site
%                 in the system through to the specified site
%
% INPUTS
% siteTensor    : complex double array, rank-4, the dmpo tensor for the
%                 specified site -- to be included in the contraction
%                 from the first site in the system
% mpo           : complex double array, rank-6, mpo tensor for the operator
%                 which is being evaluated in <rho|| O ||rho>
% leftBlock     : complex double array, rank-3, the contraction from the first
%                 site through to the site immediately before the specified site
%                 -- if the specified site is the first, then leftBlock = 1
% ROW_SIZE      : integer, the size of the first virtual dimension of siteTensor
% COL_SIZE      : integer, the size of the second virtual dimension of
%                 siteTensor
% HILBY         : integer, the size of the physical dimensions of siteTensor
% OP_COL        : integer, the size of the second virtual dimension of mpo

function [updateBlock] = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE, HILBY, OP_COL)
    updateBlock = zeros(COL_SIZE, COL_SIZE, OP_COL);

    conjTensor = zeros(ROW_SIZE, COL_SIZE, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            conjTensor(:, :, bra, ket) = conj(siteTensor(:, :, bra, ket));
        end
    end

    % permute input arrays
    % original indexing:
    %   siteTensor(row, col, bra, ket)
    %   conjTensor(conjRow, conjCol, conjBra, conjKet)
    %   mpo(conjBra, conjKet, bra, ket, opRow, opCol)
    %   leftBlock(conjCol, opRow, row)
    conjTensor = permute(conjTensor, [2, 3, 1, 4]);
    mpo = permute(mpo, [1, 5, 2, 3, 4, 6]);
    leftBlock = permute(leftBlock, [2, 3, 1]);

    for opCol = 1 : 1 : OP_COL
        AWFA = 0;
        for conjKet = 1 : 1 : HILBY
            for conjCol = 1 : 1 : ROW_SIZE
                WFA = 0;
                for bra = 1 : 1 : HILBY
                    for ket = 1 : 1 : HILBY
                        WFA = WFA + mpo(:, :, conjKet, bra, ket, opCol) * leftBlock(:, :, conjCol) * siteTensor(:, :, bra, ket);
                    end
                end
                AWFA = AWFA + conjTensor(:, :, conjCol, conjKet) * WFA;
            end
        end
        updateBlock(:, :, opCol) = AWFA;
    end

    updateBlock = permute(updateBlock, [1, 3, 2]);
end
