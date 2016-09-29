% EffLSparse.m
% function which returns the effective Liouvillian for a particular site
% in the form of a matrix, ready to be eigensolved, builds the Liouvillian
% matrix directly as a sparse matrix
% Oliver Thomson Brown
% 2016-09-28
%
% [effectiveLiouv] = EffLSparse(TARGET, dmpo, mpo, left, right)
%
% RETURN
% effectiveLiouv    : sparse complex double matrix, the effective Liovillian
%                     for the site specified by TARGET, reshaped to be a matrix
%
% INPUTS
% TARGET            : integer, the site for which the effective Liouvillian
%                     should be formed
% dmpo              : cell array, a density matrix product operator
% mpo               : cell array, mpo{n} contains the 6-dimensional matrix
%                     product operator tensor for the site n
% left              : cell array, contains the left blocks for each site
% right             : cell array, contains the right blocks for each site

function [effectiveLiouv] = EffLSparse(TARGET, dmpo, mpo, left, right)
    [ROW_SIZE, COL_SIZE, HILBY, ~] = size(dmpo{TARGET});
    LDIM = ROW_SIZE * COL_SIZE * HILBY^2;
    siteMPO = mpo{TARGET};
    lBlock = left{TARGET};
    rBlock = right{TARGET};

    % permute input arrays to avoid for-loops
    % original indexing:
    %   lBlock(conjCol, opRow, row)
    %   siteMPO(conjBra, conjKet, bra, ket, opRow, opCol)
    %   rBlock(conjRow, opCol, col)
    lBlock = permute(lBlock, [3, 2, 1]);
    siteMPO = permute(siteMPO, [5, 6, 1, 2, 3, 4]);
    rBlock = permute(rBlock, [2, 3, 1]);

    lRow = 1;
    lCol = 1;
    lVal = 0;
    for conjKet = 1 : 1 : HILBY
        for conjBra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                for bra = 1 : 1 : HILBY
                    % check if our mpo slice is non-zero
                    mpoSlice = siteMPO(:, :, conjBra, conjKet, bra, ket);
                    if any(mpoSlice(:))
                        chunkCol = (ket - 1) * HILBY * ROW_SIZE * COL_SIZE ...
                                    + (bra - 1) * ROW_SIZE * COL_SIZE;
                        for conjRow = 1 : 1 : COL_SIZE
                            for conjCol = 1 : 1 : ROW_SIZE
                                chunkRow = ...
                                (conjKet - 1) * HILBY * ROW_SIZE * COL_SIZE ...
                                + (conjBra - 1) * ROW_SIZE * COL_SIZE ...
                                + (conjCol - 1) * COL_SIZE + (conjRow - 1);
                                % put together a 2d slice
                                chunk = lBlock(:, :, conjCol) * mpoSlice ...
                                        * rBlock(:, :, conjRow);

                                % make sure the zeros are zero
                                chunk(abs(chunk) < eps) = 0;
                                chunk = reshape(transpose(chunk), ...
                                                [1, ROW_SIZE*COL_SIZE]);

                                [cRow, cCol, cVal] = find(chunk);
                                cRow = cRow + chunkRow;
                                cCol = cCol + chunkCol;
                                lRow = [lRow, cRow];
                                lCol = [lCol, cCol];
                                lVal = [lVal, cVal];
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
