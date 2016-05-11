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

	% bring to right-canonical
	compDMPO = RCan(compDMPO, LENGTH : -1 : 2);

	for site = 1 : 1 : LENGTH - 1
		compDMPO = LCan(compDMPO, site);
		[rowSz, colSz, ~, ~] = size(compDMPO{site});
		% only need to modify tensors that are too large
		if colSz > COMPRESS || rowSz > COMPRESS
			M = reshape(compDMPO{site}, [rowSz, colSz, HILBY^2]);
	        M = permute(M, [1, 3, 2]);
	        M = reshape(M, [rowSz * HILBY^2, colSz]);

			% SVD decomposition
			[U, S, V] = svd(M, 0);

			uCol = min(size(U,2), COMPRESS);
			sDim = min(size(S,1), COMPRESS);

			M = U(:, 1 : uCol) * S(1 : sDim, 1 : sDim);

			rowSz = min(size(compDMPO{site}, 1), COMPRESS);
			colSz = min(size(compDMPO{site}, 2), COMPRESS);

			% reshape back to site tensor format
			M = reshape(M, [rowSz, HILBY^2, colSz]);
			M = permute(M, [1, 3, 2]);
			compDMPO{site} = reshape(M, [rowSz, colSz, HILBY, HILBY]);

			% next site along
			V = ctranspose(V);
			vRow = min(size(V,1), COMPRESS);
			nCol = size(compDMPO{site + 1}, 2);
			V = V(1 : vRow, :);
			N = zeros(vRow, nCol, HILBY, HILBY);
			for bra = 1 : 1 : HILBY
				for ket = 1 : 1 : HILBY
					N(:, :, bra, ket) = V * compDMPO{site + 1}(:, :, bra, ket);
				end
			end
			compDMPO{site + 1} = N;
		end
	end
end
