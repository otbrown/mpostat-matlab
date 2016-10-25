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
% route     : integer array, the sites which should be brought into
%             left-canonical form, in increasing order -- note that the last
%             site cannot be the last site in the system, as the routine
%             multiplies into the next site along. If a single integer is
%             supplied, only that site is into canonical form, but the next
%             site along is still affected

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
