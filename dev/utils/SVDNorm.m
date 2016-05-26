% SVDNorm.m
% normalises the tensors in a density matrix product operator
% NOTE: This is strictly not the same as  normalising the density matrix itself
% it WILL alter the density matrix and should only be used during initialisation
% of an arbitrary state, not for normalisation
% Oliver Thomson Brown
% 2016-02-16

function [normDMPO] = SVDNorm(dmpo)
    % allocate return
    normDMPO = dmpo;

    % gather constants
    LENGTH = size(normDMPO, 1);
    HILBY = size(normDMPO{1}, 3);

    for site = 1 : 1 : LENGTH - 1
        [rowSz, colSz, ~, ~] = size(normDMPO{site});
        M = reshape(normDMPO{site}, [rowSz, colSz, HILBY^2]);
        M = permute(M, [1, 3, 2]);
        M = reshape(M, [rowSz * HILBY^2, colSz]);

        [U, S, V] = svd(M, 'econ');
        V = ctranspose(V);
        dim = size(U, 2);

        U2 = zeros(rowSz * HILBY^2, colSz);
        U2(:, 1 : dim) = U;
        U2 = reshape(U2, [rowSz, HILBY^2, colSz]);
        U2 = permute(U2, [1, 3, 2]);
        normDMPO{site} = reshape(U2, [rowSz, colSz, HILBY, HILBY]);

        rowSz = colSz;
        colSz = size(normDMPO{site + 1}, 2);

        S2 = zeros(rowSz, dim);
        S2(1 : dim, :) = S;

        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                normDMPO{site + 1}(:, :, bra, ket) = ...
                S2 * V * normDMPO{site + 1}(:, :, bra, ket);
            end
        end
    end

    [rowSz, colSz, ~, ~] = size(normDMPO{LENGTH});
    M = reshape(normDMPO{LENGTH}, [rowSz, colSz, HILBY^2]);
    M = permute(M, [1, 3, 2]);
    M = reshape(M, [rowSz * HILBY^2, colSz]);

    [U, S, V] = svd(M, 'econ');
    U2 = zeros(rowSz * HILBY^2, 1);
    U2(:, 1) = U;
    U2 = reshape(U2, [rowSz, HILBY^2, 1]);
    U2 = permute(U2, [1, 3, 2]);
    normDMPO{LENGTH} = reshape(U2, [rowSz, 1, HILBY, HILBY]);
end
