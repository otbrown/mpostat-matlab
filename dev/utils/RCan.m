% RCan.m
% Oliver Thomson Brown
% 2016-03-16
%
% rdmpo = RCan(dmpo, route)
%
% RETURN
% rdmpo     : cell array, the dmpo in right-canonical form along the specified
%             route
%
% INPUTS
% dmpo      : cell array, a density-matrix product operator
% route     : integer array, the sites which should be brought into
%             right-canonical form, in decreasing order. The last site cannot be
%             the first site in the system since the next site along is also
%             affected by this procedure

function [rdmpo] = RCan(dmpo, route)
    % gather constants
    LENGTH = size(dmpo, 1);
    HILBY = size(dmpo{1}, 3);

    if route(end) <= 1
		msgID = 'RCan:BadRoute';
		msg = sprintf('Route cannnot extend to (or exceed) the first site in the system. Your route ended at %d.', route(end));
		badRouteException = MException(msgID, msg);
		throw(badRouteException);
    end

    % allocate return
    rdmpo = dmpo;

    for site = route
        % manipulate site tensor into matrix
        [rowSz, colSz, ~, ~] = size(rdmpo{site});
        M = reshape(rdmpo{site}, [rowSz, colSz * HILBY^2]);

        % SVD decomposition
        % ISSUE! This uses a lot of mem compared to the LCan QR
        [U, S, V] = svd(M, 0);
        V = ctranspose(V);

        rowSz = size(V, 1);

        % manipulate V back into a rank-4 tensor and embed
        rdmpo{site} = reshape(V, [rowSz, colSz, HILBY, HILBY]);

        % multiply U * S into the next site along
        rowSz = size(rdmpo{site - 1}, 1);
        colSz = size(S, 2);
        N = zeros(rowSz, colSz, HILBY, HILBY);
        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                N(:, :, bra, ket) = rdmpo{site - 1}(:, :, bra, ket) * U * S;
            end
        end
        rdmpo{site - 1} = N;
    end
end
