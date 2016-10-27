% DMPOCompress.m
% truncates the virtual dimensions of a density matrix product operator to
% be no greater than COMPRESS, SVD norms in the process
% THIS IS NOT A TRACE PRESERVING OPERATION
% Oliver Thomson Brown
% 2016-02-08
%
% [ compDMPO ] = DMPOCompress(dmpo, COMPRESS)
%
% RETURN
% compDMPO:		cell array, compressed dmpo
%
% INPUT
% dmpo:			cell array, uncompressed dmpo
% COMPRESS:		integer, the maximum matrix dimension to be allowed in the
%				compressed dmpo

function [compDMPO] = DMPOCompress(dmpo, COMPRESS)
	% gather constants
	LENGTH = size(dmpo, 1);
	MIDSITE = ceil(LENGTH / 2);
	[~, OLD_COMPRESS, HILBY, ~] = size(dmpo{MIDSITE});

	if COMPRESS < HILBY^2
		msgID = 'DMPOCompress:BadCOMPRESS';
		msg = sprintf('Minimum matrix dimension is %d. Supplied COMPRESS value was %d.', HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
	elseif COMPRESS > OLD_COMPRESS
		msgID = 'DMPOCompress:BadCOMPRESS';
		msg = sprintf('Virtual dimension growth not currently supported. Current maximum virtual dimension is: %g.', OLD_COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
	end

	% allocate return
	compDMPO = dmpo;

	% bring to mixed-canonical
	compDMPO = Can(compDMPO, LENGTH : -1 : 2, 'R');
	compDMPO = Can(compDMPO, 1, 'L');

	[rowSz, colSz, ~, ~] = size(compDMPO{1});

	for site = 1 : 1 : LENGTH - 1
		% only need to modify tensors that are too large
		if colSz > COMPRESS
			siteTensor = reshape(compDMPO{site},[rowSz, colSz, HILBY^2]);
	        siteTensor = permute(siteTensor, [1, 3, 2]);
	        siteTensor = reshape(siteTensor, [rowSz * HILBY^2, colSz]);

			% SVD decomposition
			[U, S, V] = svds(siteTensor, COMPRESS);
			U = reshape(U, [rowSz, HILBY^2, COMPRESS]);
			U = permute(U, [1, 3, 2]);
			compDMPO{site} = reshape(U, [rowSz, COMPRESS, HILBY, HILBY]);

			% next site along
			rowSz = COMPRESS;
			colSz = size(compDMPO{site + 1}, 2);
			SVNextSite = zeros(rowSz, colSz, HILBY, HILBY);
			SV = S * ctranspose(V);
			for bra = 1 : 1 : HILBY
				for ket = 1 : 1 : HILBY
					SVNextSite(:, :, bra, ket) = SV  ...
						* compDMPO{site + 1}(:, :, bra, ket);
				end
			end
			compDMPO{site + 1} = SVNextSite;
		else
			compDMPO = Can(compDMPO, site, 'L');
			[rowSz, colSz, ~, ~] = size(compDMPO{site + 1});
		end
	end

	% use last site to SVD norm the compressed dmpo
	siteTensor = reshape(compDMPO{LENGTH}, [rowSz, colSz, HILBY^2]);
	siteTensor = permute(siteTensor, [1, 3, 2]);
	siteTensor = reshape(siteTensor, [rowSz * HILBY^2, colSz]);

	[U, ~, ~] = svd(siteTensor, 0);
	U = reshape(U, [rowSz, HILBY^2, 1]);
	U = permute(U, [1, 3, 2]);
	compDMPO{LENGTH} = reshape(U, [rowSz, 1, HILBY, HILBY]);
end
