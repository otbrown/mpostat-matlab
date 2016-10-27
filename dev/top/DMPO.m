% DMPO.m
% function which generates an arbitrary complex valued density matrix
% product operator, which is trace-normed
% Oliver Thomson Brown
% 2016-02-05
%
% [ dmpo ] = DMPO(HILBY, LENGTH, COMPRESS)
%
% RETURN
% dmpo:	cell, an arbitrary density matrix product operator
%
% INPUT
% HILBY:	integer, size of the local state space
% LENGTH:	integer, the number of sites in the system
% COMPRESS:	integer, the maximum virtual dimension of any given tensor

function [dmpo] = DMPO(HILBY, LENGTH, COMPRESS)

	% COMPRESS == 0 means no compression
	if COMPRESS == 0
		COMPRESS = Inf;
	elseif COMPRESS < HILBY^2
		msgID = 'DMPO:BadCOMPRESS';
		msg = sprintf(['Minimum matrix dimension is %d. Supplied ', ...
						'COMPRESS value was %d.'], HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
	end

	% allocate return
	dmpo = cell(LENGTH, 1);

	% first and last site
	Z = (1 + 1i) / sqrt(2);
	dmpo{1} = Z * ones(1, HILBY^2, HILBY, HILBY);
	dmpo{LENGTH} = Z * ones(HILBY^2, 1, HILBY, HILBY);

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

		dmpo{site} = Z * ones(rowSz, colSz, HILBY, HILBY);
	end

	% trace norm
	dmpo = TrNorm(dmpo);
end
