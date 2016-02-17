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

	% THIS GIVES CORRECT SHAPE BUT DRASTICALLY ALTERS TRACE
	for site = 1 : 1 : LENGTH - 1
		[rowSz, colSz, ~, ~] = size(compDMPO{site});

		

		M = reshape(compDMPO{site}, [rowSz, colSz, HILBY^2]);
		M = permute(M, [1, 3, 2]);
		M = reshape(M, [rowSz * HILBY^2, colSz]);

		[Q, R] = qr(M, 0);

		colSz = size(Q,2);

		Q = reshape(Q, [rowSz, HILBY^2, colSz]);
		Q = permute(Q, [1, 3, 2]);
		Q = reshape(Q, [rowSz, colSz, HILBY, HILBY]);

		% actual compression kicks in here
		rowSz = min(rowSz, COMPRESS);
		colSz = min(colSz, COMPRESS);
		compDMPO{site} = Q(1 : rowSz, 1 : colSz, :, :);

		rowSz = size(R, 1);
		colSz = size(compDMPO{site + 1}, 2);
		N = zeros(rowSz, colSz, HILBY, HILBY);
		for bra = 1 : 1 : HILBY
			for ket = 1 : 1 : HILBY
				N(:, :, bra, ket) = R * compDMPO{site + 1}(:, :, bra, ket);
			end
		end
		compDMPO{site + 1} = N;
	end
end
