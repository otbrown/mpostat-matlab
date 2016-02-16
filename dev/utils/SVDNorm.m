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

        [Q, R] = qr(M, 0);

        Q = reshape(Q, [rowSz, HILBY^2, colSz]);
        Q = permute(Q, [1, 3, 2]);
        normDMPO{site} = reshape(Q, [rowSz, colSz, HILBY, HILBY]);

        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                normDMPO{site + 1}(:, :, bra, ket) = R * normDMPO{site + 1}(:, :, bra, ket);
            end
        end
    end

    [rowSz, colSz, ~, ~] = size(normDMPO{LENGTH});
    M = reshape(normDMPO{LENGTH}, [rowSz, colSz, HILBY^2]);
    M = permute(M, [1, 3, 2]);
    M = reshape(M, [rowSz * HILBY^2, colSz]);

    [Q, ~] = qr(M, 0);

    Q = reshape(Q, [rowSz, HILBY^2, colSz]);
    Q = permute(Q, [1, 3, 2]);
    normDMPO{LENGTH} = reshape(Q, [rowSz, colSz, HILBY, HILBY]);

end
