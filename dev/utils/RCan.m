% RCan.m
% function which returns a site tensor in right-canonical form and its
% neighbour which has the norm carried into it
% Oliver Thomson Brown
% 2016-03-16
%
% [ canSite, nextSiteUS ] = RCan(siteTensor, nextSiteTensor, HILBY, ...
%                                   ROW_SIZE, COL_SIZE, NEXT_ROW)
%
% RETURN
% canSite:      complex double, rank-4, the site tensor which has been made
%               right-canonical
% nextSiteUS:   complex double, rank-4, the site tensor which follows
%               canSite in the system -- has been modified in the process
%               of making canSite right-canonical
%
% INPUT
% siteTensor:       complex double, rank-4, the site tensor which is to be
%                   made right-canonical
% nextSiteTensor:   complex double, rank-4, the site tensor for the site
%                   which follows siteTensor
% HILBY:            integer, the physical dimension of the system
% ROW_SIZE:         integer, the first virtual dimension of siteTensor
% COL_SIZE:         integer, the second virtual dimension of siteTensor
% NEXT_ROW:         integer, the first virtual dimension of nextSiteTensor

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
