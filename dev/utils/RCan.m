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
%             right-canonical form, in decreasing order. The last site cannot
%             be the first site in the system since the next site along is
%             also affected by this procedure

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

    [rowSz, colSz, ~, ~] = size(rdmpo{route(1)});

    for site = route
        % manipulate site tensor into matrix
        M = reshape(rdmpo{site}, [rowSz, colSz * HILBY^2]);

        % SVD decomposition
        [U, S, V] = svd(M, 'econ');
        vr = size(V, 2);
        V2 = zeros(rowSz, colSz * HILBY^2);
        V2(1 : vr, :) = ctranspose(V);

        % manipulate V back into a rank-4 tensor and embed
        rdmpo{site} = reshape(V2, [rowSz, colSz, HILBY, HILBY]);

        % multiply U * S into the next site along
        colSz = rowSz;
        rowSz = size(rdmpo{site - 1}, 1);

        US = zeros(colSz);
        US(:, 1 : vr) = U * S;

        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                rdmpo{site - 1}(:, :, bra, ket) = ...
                rdmpo{site - 1}(:, :, bra, ket) * US;
            end
        end
    end
end
