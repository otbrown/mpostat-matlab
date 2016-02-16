% DMPO.m
% function which generates an arbitrary density matrix product operator
% Oliver Thomson Brown
% 2016-02-05

function [dmpo] = DMPO(HILBY, LENGTH, COMPRESS)

	% COMPRESS == 0 means no compression
	if COMPRESS == 0
		COMPRESS = Inf;
	end

	% allocate return
	dmpo = cell(LENGTH, 1);

	% first and last site
	dmpo{1} = rand(1, HILBY^2, HILBY, HILBY) + 1i*rand(1, HILBY^2, HILBY, HILBY);
	dmpo{LENGTH} = rand(HILBY^2, 1, HILBY, HILBY) + 1i*rand(HILBY^2, 1, HILBY, HILBY);

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

		dmpo{site} = rand(rowSz, colSz, HILBY, HILBY) + 1i*rand(rowSz, colSz, HILBY, HILBY);
	end

	% SVD normalise -- This normalises the tensors NOT THE STATE
	dmpo = SVDNorm(dmpo);

	% make Hermitian
	% *** BUILD THIS ***
	% rho = (rho + rho^H)/2
	% dmpo = Herm(dmpo);

	% trace norm
	% *** BUILD THIS ***
	% dmpo = TrNorm(dmpo);
end
