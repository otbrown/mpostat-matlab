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
    updateBlock = zeros(ROW_SIZE, ROW_SIZE, OP_ROW);

    conjTensor = zeros(ROW_SIZE, COL_SIZE, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            conjTensor(:, :, bra, ket) = conj(siteTensor(:, :, bra, ket));
        end
    end

    % permute input arrays to aid optimisation
    % original indexing:
    %   siteTensor(row, col, bra, ket)
    %   rightBlock(conjRow, opCol, col)
    %   mpo(conjBra, conjKet, bra, ket, opRow, opCol)
    %   conjTensor(conjRow, conjCol, conjBra, conjKet)
    siteTensor = permute(siteTensor, [2, 1, 3, 4]);
    rightBlock = permute(rightBlock, [2, 3, 1]);
    mpo = permute(mpo, [1, 6, 2, 3, 4, 5]);
    conjTensor = permute(conjTensor, [1, 3, 2, 4]);

    for opRow = 1 : 1 : OP_ROW
        BWFB = 0;
        for conjKet = 1 : 1 : HILBY
            for conjRow = 1 : 1 : COL_SIZE
                WFB = 0;
                for bra = 1 : 1 : HILBY
                    for ket = 1 : 1 : HILBY
                        WFB = WFB + mpo(:, :, conjKet, bra, ket, opRow) * rightBlock(:, :, conjRow) * siteTensor(:, :, bra, ket);
                    end
                end
                BWFB = BWFB + conjTensor(:, :, conjRow, conjKet) * WFB;
            end
        end
        updateBlock(:, :, opRow) = BWFB;
    end

    updateBlock = permute(updateBlock, [1, 3, 2]);
end
