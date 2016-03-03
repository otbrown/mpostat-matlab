% FWBase.m
% function to return an array containing the fixed width
% base b representation of a number, n
% Oliver Thomson Brown
% 2016-02-03

function [bits] = FWBase(n, BASE, WIDTH)
    % check that WIDTH is high enough
    minWidth = floor(log(n) / log(BASE)) + 1;
    if n < 0 || mod(n,1) ~= 0
        msgID = 'FWBase:BadN';
        msg = sprintf('FWBase currently only handles positive integers.');
        nException = MException(msgID, msg);

        throw(nException);

    elseif WIDTH < minWidth
        msgID = 'FWBase:BadWidth';
        msg = sprintf('%d needs at least %d bits in base %d. WIDTH: %d.', n, minWidth, BASE, WIDTH);
        widthException = MException(msgID, msg);

        throw(widthException);
    end

    bits = zeros(WIDTH,1);
    for i = WIDTH : -1 : (WIDTH - minWidth + 1)
        bits(i) = mod(n, BASE);
        n = fix(n / BASE);
    end
end
