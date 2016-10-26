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

function [canSite, nextSiteUS] = RCan(siteTensor, nextSiteTensor, HILBY, ROW_SIZE, COL_SIZE, NEXT_ROW)
    siteTensor = reshape(siteTensor, [ROW_SIZE, COL_SIZE * HILBY^2]);

    % SVD decomposition
    [U, S, V] = svd(siteTensor, 'econ');

    % manipulate V back into a rank-4 tensor and embed
    canSite = reshape(ctranspose(V), [ROW_SIZE, COL_SIZE, HILBY, HILBY]);

    % multiply U * S into the next site along
    US = U * S;

    nextSiteUS = zeros(NEXT_ROW, ROW_SIZE, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            nextSiteUS(:, :, bra, ket) = ...
            nextSiteTensor(:, :, bra, ket) * US;
        end
    end
end
