% GrowLeft.m
% contracts a given dmpo site tensor into the left block -- a contraction from the
% first site through to the site prior to this one.
% NAIEVE IMPLEMENTATION -- OPTIMISE!
% Oliver Thomson Brown
% 2016-03-22
%
% updateBlock = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
%
% RETURN
% updateBlock   : complex double array, rank-3, contraction from the first site in %                 the system through to the specified site
%
% INPUTS
% siteTensor    : complex double array, rank-4, the dmpo tensor for the specified
%                 site -- to be included in the contraction from the first site in %                 the system
% mpo           : complex double array, rank-6, mpo tensor for the operator which
%                 is being evaluated in <rho|| O ||rho>
% leftBlock     : complex double array, rank-3, the contraction from the first
%                 site through to the site immediately before the specified site
%                 -- if the specified site is the first, then leftBlock = 1
% ROW_SIZE      : integer, the size of the first virtual dimension of siteTensor
% COL_SIZE      : integer, the size of the second virtual dimension of siteTensor
% HILBY         : integer, the size of the physical dimensions of siteTensor
% OP_ROW        : integer, the size of the first virtual dimension of mpo
% OP_COL        : integer, the size of the second virtual dimension of mpo

function [updateBlock] = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)

updateBlock = zeros(COL_SIZE, OP_COL, COL_SIZE);
conjTensor = zeros(COL_SIZE, ROW_SIZE, HILBY, HILBY);
for bra = 1 : 1 : HILBY
    for ket = 1 : 1 : HILBY
        conjTensor(:, :, bra, ket) = ctranspose(siteTensor(:, :, bra, ket));
    end
end

for conjRow = 1 : 1 : COL_SIZE
    for opCol = 1 : 1 : OP_COL
        for col = 1 : 1 : COL_SIZE
            AWFA = 0;
            for conjBra = 1 : 1 : HILBY
                for conjKet = 1 : 1 : HILBY
                    for conjCol = 1 : 1 : ROW_SIZE
                        WFA = 0;
                        for bra = 1 : 1 : HILBY
                            for ket = 1 : 1 : HILBY
                                for opRow = 1 : 1 : OP_ROW
                                    FA = 0;
                                    for row = 1 : 1 : ROW_SIZE
                                        FA = FA + leftBlock(conjCol, opRow, row) * siteTensor(row, col, bra, ket);
                                    end
                                    WFA = WFA + mpo(bra, ket, conjBra, conjKet, opRow, opCol) * FA;
                                end
                            end
                        end
                        AWFA = AWFA + conjTensor(conjRow, conjCol, conjBra, conjKet) * WFA;
                    end
                end
            end
            updateBlock(conjRow, opCol, col) = AWFA;
        end
    end
end

end
