% DMPOEnlarge.m
% increases the virtual dimensions of an existing dmpo
% the extra space is padded with zeros
% Oliver Thomson Brown
% 2016-11-16
%
% [ bigDMPO ] = DMPOEnlarge(dmpo, COMPRESS, HILBY, LENGTH)
%
% RETURN
% bigDMPO:  cell, the density matrix product operator with enlarged virtual
%           dimensions
%
% INPUT
% dmpo:     cell, an arbitrary density matrix product operator
% COMPRESS: integer, new maximum virtual dimension of the tensors of dmpo
% HILBY:    integer, the size of the local state space
% LENGTH:   integer, the number of sites in the system

function [bigDMPO] = DMPOEnlarge(dmpo, COMPRESS, HILBY, LENGTH)
    bigDMPO = dmpo;

    rowSz = 1;
    for site = 1 : 1 : LENGTH
        if site < ceil(LENGTH / 2);
            len = site;
        else
            len = LENGTH - site;
        end
        colSz = min(HILBY^(2*len), COMPRESS);

        dim = size(bigDMPO{site});
        if any(dim ~= [rowSz, colSz, HILBY, HILBY])
            A = zeros(rowSz, colSz, HILBY, HILBY);
            A(1 : dim(1), 1 : dim(2), :, :) = bigDMPO{site};
            bigDMPO{site} = A;
        end

        rowSz = colSz;
    end
end
