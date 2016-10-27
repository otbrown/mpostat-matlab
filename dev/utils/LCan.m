% LCan.m
% function which returns a site tensor in left-canonical form and its
% neighbour which has the norm carried into it
% Oliver Thomson Brown
% 2016-03-15
%
% [ canSite, SVNextSite ] = LCan(siteTensor, nextSiteTensor, HILBY, ...
%                                   ROW_SIZE, COL_SIZE, NEXT_COL)
%
% RETURN
% canSite:      complex double, rank-4, the site tensor which has been
%               made left-canonical
% SVNextSite:   complex double, rank-4, the site tensor for the site which
%               follows canSite -- has been modified in the process of
%               making canSite left-canonical
%
% INPUT
% siteTensor:       complex double, rank-4, the site tensor which is to be
%                   made left-canonical
% nextSiteTensor:   complex double, rank-4, the site tensor for the site
%                   which follows siteTensor
% HILBY:            integer, physical dimension of the system
% ROW_SIZE:         integer, the first virtual dimension of siteTensor
% COL_SIZE:         integer, the second virtual dimension of siteTensor
% NEXT_COL:         integer, the second virtual dimension of nextSiteTensor

function [canSite, SVNextSite] = LCan(siteTensor, nextSiteTensor, HILBY, ROW_SIZE, COL_SIZE, NEXT_COL)
    % manipulate site tensor into a matrix
    siteTensor = reshape(siteTensor, [ROW_SIZE, COL_SIZE, HILBY^2]);
    siteTensor = permute(siteTensor, [1, 3, 2]);
    siteTensor = reshape(siteTensor, [ROW_SIZE * HILBY^2, COL_SIZE]);

    % SVD Decomposition
    [U, S, V] = svd(siteTensor, 0);

    % manipulate U into rank-4 tensor and embed in site
    U = reshape(U, [ROW_SIZE, HILBY^2, COL_SIZE]);
    U = permute(U, [1, 3, 2]);
    canSite = reshape(U, [ROW_SIZE, COL_SIZE, HILBY, HILBY]);

    % multiply SV into the next site
    SV = S * ctranspose(V);

    SVNextSite = zeros(COL_SIZE, NEXT_COL, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            SVNextSite(:, :, bra, ket) = ...
            SV * nextSiteTensor(:, :, bra, ket);
        end
    end
end
