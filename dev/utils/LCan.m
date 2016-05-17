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
    % gather constants
    LENGTH = size(dmpo, 1);
    HILBY = size(dmpo{1}, 3);

    if route(end) >= LENGTH
		msgID = 'LCan:BadRoute';
		msg = sprintf('Route cannnot extend to (or exceed) the last site in the system. System has %d sites, your route ended at %d.', LENGTH, route(end));
		badRouteException = MException(msgID, msg);
		throw(badRouteException);
    end

    % allocate return
    ldmpo = dmpo;

    for site = route
        % manipulate site tensor into a matrix
        [rowSz, colSz, ~, ~] = size(ldmpo{site});
        M = reshape(ldmpo{site}, [rowSz, colSz, HILBY^2]);
        M = permute(M, [1, 3, 2]);
        M = reshape(M, [rowSz * HILBY^2, colSz]);

        % QR decomposition
        %[Q, R] = qr(M,0);
        % SVD Decomposition
        [U, S, V] = svd(M, 'econ');
        V = ctranspose(V);
        dim = size(U, 2);

        %colSz = size(Q, 2);
        %colSz = size(US, 2);

        % manipulate Q into rank-4 tensor and embed in site
        %Q = reshape(Q, [rowSz, HILBY^2, colSz]);
        %Q = permute(Q, [1, 3, 2]);
        %ldmpo{site} = reshape(Q, [rowSz, colSz, HILBY, HILBY]);
        % reshape M to be the new site tensor
        %M = reshape(US, [rowSz, HILBY^2, colSz]);
        %M = permute(US, [1, 3, 2]);
        %ldmpo{site} = reshape(US, [rowSz, colSz, HILBY, HILBY]);
        U2 = zeros(rowSz * HILBY^2, colSz);
        U2(:, 1 : dim) = U;
        U2 = reshape(U2, [rowSz, HILBY^2, colSz]);
        U2 = permute(U2, [1, 3, 2]);
        ldmpo{site} = reshape(U2, [rowSz, colSz, HILBY, HILBY]);

        % multiply R into the next site along
        %rowSz = size(R, 1);
        %rowSz = size(V, 1);
        %colSz = size(ldmpo{site+1}, 2);
        %N = zeros(rowSz, colSz, HILBY, HILBY);
        rowSz = colSz;
        colSz = size(ldmpo{site + 1}, 2);

        S2 = zeros(rowSz, dim);
        S2(1 : dim, :) = S;

        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                %N(:, :, bra, ket) = R * ldmpo{site + 1}(:, :, bra, ket);
                %N(:, :, bra, ket) = V * ldmpo{site + 1}(:, :, bra, ket);
                ldmpo{site + 1}(:, :, bra, ket) = ...
                S2 * V * ldmpo{site + 1}(:, :, bra, ket);
            end
        end
        %ldmpo{site + 1} = N;
    end
end
