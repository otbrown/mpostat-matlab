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
	
	for site = 1 : 1 : LENGTH
		[rowSz, colSz, ~, ~] = size(dmpo{site});
		if rowSz > COMPRESS || colSz > COMPRESS
			rowSz = min(rowSz, COMPRESS);
			colSz = min(colSz, COMPRESS);
			compDMPO{site} = zeros(rowSz, colSz, HILBY, HILBY);
			% there is probably a better way than this naieve approach
			for ket = 1 : 1 : HILBY
				for bra = 1 : 1 : HILBY
					[U, S, V] = svd(dmpo{site}(:, :, bra, ket));
					U = U(1 : rowSz, 1 : rowSz);
					S = S(1 : rowSz, 1 : colSz);
					V = ctranspose(V);
					V = V(1 : colSz, 1 : colSz);
					compDMPO{site}(:, :, bra, ket) = U * S * V; 
				end
			end
		end
	end
end
