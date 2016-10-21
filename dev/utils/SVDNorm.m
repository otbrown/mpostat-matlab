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

    [rowSz, colSz, ~, ~] = size(normDMPO{1});

    for site = 1 : 1 : LENGTH - 1
        M = reshape(normDMPO{site}, [rowSz, colSz, HILBY^2]);
        M = permute(M, [1, 3, 2]);
        M = reshape(M, [rowSz * HILBY^2, colSz]);

        [U, S, V] = svd(M, 0);

        U = reshape(U, [rowSz, HILBY^2, colSz]);
        U = permute(U, [1, 3, 2]);
        normDMPO{site} = reshape(U, [rowSz, colSz, HILBY, HILBY]);

        rowSz = colSz;
        colSz = size(normDMPO{site + 1}, 2);

        SV = S * ctranspose(V);

        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                normDMPO{site + 1}(:, :, bra, ket) = ...
                SV * normDMPO{site + 1}(:, :, bra, ket);
            end
        end
    end

    M = reshape(normDMPO{LENGTH}, [rowSz, colSz, HILBY^2]);
    M = permute(M, [1, 3, 2]);
    M = reshape(M, [rowSz * HILBY^2, colSz]);

    [U, S, V] = svd(M, 0);

    U = reshape(U, [rowSz, HILBY^2, 1]);
    U = permute(U, [1, 3, 2]);
    normDMPO{LENGTH} = reshape(U, [rowSz, 1, HILBY, HILBY]);
end
