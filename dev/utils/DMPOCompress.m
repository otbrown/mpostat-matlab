% DMPOCompress.m
% truncates the virtual dimensions of a tensor network to
% be no greater than COMPRESS
% Oliver Thomson Brown
% 2016-02-08
%
% compDMPO = DMPOCompress(dmpo, COMPRESS)
%
% RETURN
% compDMPO:		cell array, compressed dmpo representation
%
% INPUTS
% dmpo:			cell array, uncompressed dmpo representation
% COMPRESS:		integer, the maximum matrix dimension to be allowed in the
%				compressed dmpo

function [compDMPO] = DMPOCompress(dmpo, COMPRESS)
	% gather constants
	LENGTH = size(dmpo, 1);
	HILBY = size(dmpo{1}, 3);

	if COMPRESS < HILBY^2
		msgID = 'DMPOCompress:BadCOMPRESS';
		msg = sprintf('Minimum matrix dimension is %d. Supplied COMPRESS value was %d.', HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
	end

	% allocate return
	compDMPO = dmpo;

	% COMPRESS COLUMN DIMENSIONS
	for site = 1 : 1 : LENGTH - 1
		[rowSz, colSz, ~, ~] = size(dmpo{site});
		if colSz > COMPRESS || rowSz > COMPRESS
			% there is probably a better way than this naieve approach
			% LIKE FOR EXAMPLE SOMETHING THAT'S NOT WRONG
			M = reshape(compDMPO{site}, [rowSz, colSz, HILBY^2]);
			M = permute(M, [1, 3, 2]);
			M = reshape(M, [rowSz*HILBY^2, colSz]);
			[U, S, V] = svd(M);
			M = U * S;
			M = reshape(M, [rowSz, colSz, HILBY^2]);
			M = permute(M, [1, 3, 2]);
			M = reshape(M, [rowSz, colSz, HILBY, HILBY]);

			rowSz = min(rowSz, COMPRESS);
			colSz = min(colSz, COMPRESS);
			compDMPO{site} = M(1 : rowSz, 1 : colSz, :, :);

			V = ctranspose(V);
			for bra = 1 : 1 : HILBY
				for ket = 1 : 1 : HILBY
					compDMPO{site + 1}(:, :, bra, ket) = V * compDMPO{site + 1}(:, :, bra, ket);
				end
			end
		end
	end

end
