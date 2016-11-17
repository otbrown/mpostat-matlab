% DMPOResize.m
% Oliver Thomson Brown
% 2016-11-17

function [rsDMPO] = DMPOResize(dmpo, COMPRESS)
    LENGTH = length(dmpo);
    MIDSITE = floor(LENGTH / 2);
    [~, OLD_COMPRESS, HILBY, ~] = size(dmpo{MIDSITE});
    EXACT = HILBY^(2*MIDSITE);

	if COMPRESS < HILBY^2
		msgID = 'DMPOResize:BadCOMPRESS';
		msg = sprintf('Minimum matrix dimension is %d. Supplied COMPRESS value was %d.', HILBY^2, COMPRESS);
		badCOMPRESSException = MException(msgID, msg);
		throw(badCOMPRESSException);
    end

    if COMPRESS < OLD_COMPRESS
        rsDMPO = DMPOCompress(dmpo, COMPRESS, HILBY, LENGTH);
        rsDMPO = TrNorm(rsDMPO);
    elseif COMPRESS > OLD_COMPRESS && OLD_COMPRESS ~= EXACT
        rsDMPO = DMPOEnlarge(dmpo, COMPRESS, HILBY, LENGTH);
    else
        rsDMPO = dmpo;
    end
end
