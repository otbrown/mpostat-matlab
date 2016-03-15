% LCan.m
% function which returns the a dmpo with the requested sites in left
% canonical form
% Oliver Thomson Brown
% 2016-03-15
%
% ldmpo = LCan(dmpo, route)
%
% RETURN
% ldmpo     : cell array, the dmpo in left-canonical form along the specified
%             route
%
% INPUTS
% dmpo      : cell array, a density matrix product operator
% route     : integer array, the sites which should be brought into left-canonical
%             form, in increasing order -- note that the last site cannot be the
%             last site in the system, as the routine multiplies into the next
%             site along. If a single integer is supplied, only that site is
%             into canonical form, but the next site along is still affected

function [ldmpo] = LCan(dmpo, route)
    % allocate return
    ldmpo = dmpo;

    % gather constants
    LENGTH = size(ldmpo, 1);
    HILBY = size(ldmpo{1}, 3);

    for site = route
        % manipulate site tensor into a matrix
        [rowSz, colSz, ~, ~] = size(ldmpo{site});
        M = reshape(ldmpo{site}, [rowSz, colSz, HILBY^2]);
        M = permute(M, [1, 3, 2]);
        M = reshape(M, [rowSz * HILBY^2, colSz]);

        % QR decomposition
        [Q, R] = qr(M,0);

        colSz = size(Q, 2);

        % manipulate Q into rank-4 tensor and embed in site
        Q = reshape(Q, [rowSz, HILBY^2, colSz]);
        Q = permute(Q, [1, 3, 2]);
        ldmpo{site} = reshape(Q, [rowSz, colSz, HILBY, HILBY]);

        % multiply R into the next site along
        rowSz = size(R, 1);
        colSz = size(ldmpo{site+1}, 2);
        N = zeros(rowSz, colSz, HILBY, HILBY);
        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                N(:, :, bra, ket) = R * ldmpo{site + 1}(:, :, bra, ket);
            end
        end
        ldmpo{site + 1} = N;
    end
end
