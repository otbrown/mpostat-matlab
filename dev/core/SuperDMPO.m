% SuperDMPO.m
% (Super)position (D)ensity (M)atrix (P)roduct (O)perator
% forms a matrix product density operator with every tensor as a single element
% identity (a one in the top left), which is then trace-normed
% represents an even mix of every possible state
% Oliver Thomson Brown
% 2016-03-14
%
% dmpo = SuperDMPO(HILBY, LENGTH, COMPRESS)
%
% RETURN
% dmpo      : cell array, the density matrix product operator
%
% INPUTS
% HILBY		: integer, size of the local state space
% LENGTH	: integer, the number of sites in the system
% COMPRESS	: integer, the maximum virtual dimension of any given tensor

function [dmpo] = SuperDMPO(HILBY, LENGTH, COMPRESS)
    % COMPRESS == 0 means no compression
    if COMPRESS == 0
        COMPRESS = Inf;
    elseif COMPRESS < HILBY^2
		msgID = 'MDMPO:BadCOMPRESS';
		msg = sprintf('Minimum matrix dimension is %d. Supplied COMPRESS value was %d.', HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
	end

    % allocate return
    dmpo = cell(LENGTH, 1);

    % first and last site
    dmpo{1} = zeros(1, HILBY^2, HILBY, HILBY);
    dmpo{1}(1, 1, :, :) = ones(1, 1, HILBY, HILBY);
    dmpo{LENGTH} = zeros(HILBY^2, 1, HILBY, HILBY);
    dmpo{LENGTH}(1, 1, :, :) = ones(1, 1, HILBY, HILBY);

    % and the rest
	colSz = HILBY^2;
	for site = 2 : 1 : LENGTH - 1
		if site < ceil(LENGTH / 2)
			len = site;
		else
			len = LENGTH - site;
		end

		rowSz = colSz;
		colSz = min(HILBY^(2*len), COMPRESS);

        dmpo{site} = zeros(rowSz, colSz, HILBY, HILBY);
        dmpo{site}(1, 1, :, :) = ones(1, 1, HILBY, HILBY);
    end

    % already definitely is Hermitian (every element of the density matrix is a 1)
    % so no need to make it Hermitian, likewise no need to SVD normalise
    % technically we should use TrNorm, but we already know fine well that the
    % trace is HILBY^LENGTH, so we can skip straight to the scalar division
    dmpo = DMPOScalarDiv(dmpo, HILBY^LENGTH);
end
