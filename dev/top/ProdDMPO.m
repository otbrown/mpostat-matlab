% ProdDMPO.m
% function which builds a dmpo representing the density matrix of a
% product state
% Oliver Thomson Brown
% 2016-02-26
%
% [ prodDMPO ] = ProdDMPO(HILBY, LENGTH, COMPRESS, STATE)
%
% RETURN
% prodDMPO: cell, a density matrix product operator representing a product
%           state of the system
%
% INPUT
% HILBY:    integer, the dimension of the local state space
% LENGTH:   integer, the size of the 1-D system
% COMPRESS: integer, the maximum dimension of the matrices, enter 0 for an
%           uncompressed MPS
% STATE:    integer, the decimal representation of the product state you
%           wish to form

function [prodDMPO] = ProdDMPO(HILBY, LENGTH, COMPRESS, STATE)
    % COMPRESS == 0 means no compression
    if COMPRESS == 0
        COMPRESS = Inf;
    elseif COMPRESS < HILBY^2
		msgID = 'ProdDMPO:BadCOMPRESS';
		msg = sprintf(['Minimum matrix dimension is %d. Supplied ', ...
                        'COMPRESS value was %d.'], HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
    end

    % generate indices for state
    local = FWBase(STATE, HILBY, LENGTH) + 1;

    % return allocation
    prodDMPO = cell(LENGTH, 1);

    % first and last site
    prodDMPO{1} = zeros(1, HILBY^2, HILBY, HILBY);
    prodDMPO{1}(:, :, local(1), local(1)) = eye(1, HILBY^2);
    prodDMPO{LENGTH} = zeros(HILBY^2, 1, HILBY, HILBY);
    prodDMPO{LENGTH}(:, :, local(LENGTH), local(LENGTH)) = eye(HILBY^2, 1);

    % all the other sites
    colSz = HILBY^2;

    for site  = 2 : 1 : LENGTH - 1
        if site < ceil(LENGTH / 2)
            len = site;
        else
            len = LENGTH - site;
        end

        rowSz = colSz;
        colSz = min(HILBY^(2*len), COMPRESS);

        prodDMPO{site} = zeros(rowSz, colSz, HILBY, HILBY);
        prodDMPO{site}(:, :, local(site), local(site)) = eye(rowSz, colSz);
    end
end
