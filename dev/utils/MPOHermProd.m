% MPOHermProd.m
% function which returns the MPO representation of ctranspose(O)*O when
% supplied with the mpo representation of O
% Oliver Thomson Brown
% 2016-09-30
%
% [ hmpo ] = MPOHermProd(mpo)
%
% RETURN
% hmpo: cell, mpo representation of the product of the hermitian conjugate
%       of the supplied mpo and itself -- virtual dimensions of the mpo
%       are squared during this process
%
% INPUT
% mpo: cell, a matrix product operator representing a Liouvillian

function hmpo = MPOHermProd(mpo)
    % gather sizes
    LENGTH = length(mpo);
    [HILBY, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{2});

    % find hermitian conjugate of mpo
    dagger = cell(LENGTH, 1);
    for site = 1 : 1 : LENGTH
        dagger{site} = permute(conj(mpo{site}), [3,4,1,2,5,6]);
    end

    % build hermitian product
    hmpo = cell(LENGTH, 1);

    hmpo{1} = zeros(HILBY, HILBY, HILBY, HILBY, 1, OP_COL^2);
    for opCol1 = 1 : 1 : OP_COL
        for opCol2 = 1 : 1 : OP_COL
            jOpCol = (opCol1-1) * OP_COL + opCol2;

            A = reshape(dagger{1}(:,:,:,:,1,opCol1), [HILBY^2, HILBY^2]);
            B = reshape(mpo{1}(:,:,:,:,1,opCol2), [HILBY^2, HILBY^2]);

            hmpo{1}(:, :, :, :, 1, jOpCol) = reshape(A*B, ...
                                             [HILBY, HILBY, HILBY, HILBY]);
        end
    end

    for site = 2 : 1 : LENGTH - 1
        hmpo{site} = zeros(HILBY, HILBY, HILBY, HILBY, OP_ROW^2, OP_COL^2);
        for opRow1 = 1 : 1 : OP_ROW
            for opRow2 = 1 : 1 : OP_ROW
                jOpRow = (opRow1-1) * OP_ROW + opRow2;
                for opCol1 = 1 : 1 : OP_COL
                    for opCol2 = 1 : 1 : OP_COL
                        jOpCol = (opCol1-1) * OP_COL + opCol2;

                        A = reshape( ...
                                dagger{site}(:,:,:,:,opRow1,opCol1), ...
                                [HILBY^2, HILBY^2]);
                        B = reshape(mpo{site}(:,:,:,:,opRow2,opCol2), ...
                                [HILBY^2, HILBY^2]);

                        hmpo{site}(:,:,:,:,jOpRow,jOpCol) = ...
                        reshape(A*B, [HILBY, HILBY, HILBY, HILBY]);
                    end
                end
            end
        end
    end

    hmpo{LENGTH} = zeros(HILBY, HILBY, HILBY, HILBY, OP_ROW^2, 1);
    for opRow1 = 1 : 1 : OP_ROW
        for opRow2 = 1 : 1 : OP_ROW
            jOpRow = (opRow1-1) * OP_ROW + opRow2;

            A = reshape(dagger{LENGTH}(:,:,:,:,opRow1),[HILBY^2,HILBY^2]);
            B = reshape(mpo{LENGTH}(:,:,:,:,opRow2), [HILBY^2, HILBY^2]);

            hmpo{LENGTH}(:, :, :, :, jOpRow) = ...
            reshape(A*B, [HILBY, HILBY, HILBY, HILBY]);
        end
    end
end
