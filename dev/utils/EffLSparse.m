% EffLSparse.m
% function which returns the effective Liouvillian for a particular site
% in the form of a matrix, ready to be eigensolved, builds the Liouvillian
% matrix directly as a sparse matrix
% Oliver Thomson Brown
% 2016-09-28
%
% [ effectiveLiouv ] = EffLSparse(lBlock, siteMPO, rBlock, ROW_SIZE, ...
%                                   COL_SIZE, HILBY)
%
% RETURN
% effectiveLiouv:   sparse complex double, the effective Liouvillian for
%                   the site specified by TARGET
%
% INPUT
% lBlock:   complex double, rank-3 left block tensor which contains the
%           contraction through the system from the first site up to the
%           target site
% siteMPO:  complex double, rank-6 matrix product operator tensor for the
%           target site
% rBlock:   complex double, rank-3 right block tensor which contains the
%           contraction through the sytem from the target site to the last
% ROW_SIZE: double, the size of the first virtual dimension of the density
%           matrix product operator for the target site
% COL_SIZE: double, the size of the second virtual dimension of the
%           density matrix product operator for the target site
% HILBY:    double, the size of the local state space

function [effectiveLiouv] = EffLSparse(lBlock, siteMPO, rBlock, ROW_SIZE, COL_SIZE, HILBY)
    LDIM = ROW_SIZE * COL_SIZE * HILBY^2;

    % permute input arrays to avoid for-loops
    % original indexing:
    %   lBlock(conjCol, opRow, row)
    %   siteMPO(conjBra, conjKet, bra, ket, opRow, opCol)
    %   rBlock(conjRow, opCol, col)
    lBlock = permute(lBlock, [3, 2, 1]);
    siteMPO = permute(siteMPO, [5, 6, 1, 2, 3, 4]);
    rBlock = permute(rBlock, [2, 3, 1]);

    numNZ = ceil(0.1 * LDIM^2);
    growNZ = ceil(0.05 * LDIM^2);
    lRow = ones(numNZ, 1);
    lCol = ones(numNZ, 1);
    lVal = zeros(numNZ, 1);
    ldex = 0;

    for conjKet = 1 : 1 : HILBY
        for conjBra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                for bra = 1 : 1 : HILBY
                    % check if our mpo slice is non-zero
                    mpoSlice = siteMPO(:, :, conjBra, conjKet, bra, ket);
                    if any(mpoSlice(:))
                        chunkCol = (ket - 1) * HILBY * ROW_SIZE ...
                                    * COL_SIZE + (bra - 1) * ROW_SIZE ...
                                    * COL_SIZE;
                        for conjRow = 1 : 1 : COL_SIZE
                            for conjCol = 1 : 1 : ROW_SIZE
                                chunkRow = ...
                                (conjKet - 1) * HILBY * ROW_SIZE * ...
                                COL_SIZE + (conjBra - 1) * ROW_SIZE * ...
                                COL_SIZE + (conjCol - 1) * COL_SIZE + ...
                                (conjRow - 1);

                                % put together a 2d slice
                                chunk = lBlock(:, :, conjCol) ...
                                        * mpoSlice * rBlock(:, :, conjRow);
                                chunk = reshape(transpose(chunk), ...
                                                [1, ROW_SIZE*COL_SIZE]);

                                chunk(abs(chunk) < eps) = 0;
                                [cRow, cCol, cVal] = find(chunk);
                                cRow = cRow + chunkRow;
                                cCol = cCol + chunkCol;

                                chunkSize = length(cVal);
                                if (ldex + chunkSize) > numNZ
                                    lRow = [lRow; ones(growNZ, 1)];
                                    lCol = [lCol; ones(growNZ, 1)];
                                    lVal = [lVal; zeros(growNZ, 1)];
                                    numNZ = numNZ + growNZ;
                                end

                                lRow(ldex + 1 : ldex + chunkSize) = cRow;
                                lCol(ldex + 1 : ldex + chunkSize) = cCol;
                                lVal(ldex + 1 : ldex + chunkSize) = cVal;

                                ldex = ldex + chunkSize;
                            end
                        end
                    end
                end
            end
        end
    end

    % build and return
    effectiveLiouv = sparse(lRow, lCol, lVal, LDIM, LDIM);
end
