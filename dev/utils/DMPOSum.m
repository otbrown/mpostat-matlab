% DMPOSum.m
% computes the sum of to density matrix product operators
% Oliver Thomson Brown
% 2016-02-08
%
% sumDMPO = DMPOSum(rhoA, rhoB)
%
% RETURN
% sumDMPO:			cell array, rhoA + rhoB
%
% INPUTS
% rhoA:			cell array, a density matrix expressed as an MPO
% rhoB:			cell array, a density matrix expressed as an MPO

function [sumDMPO] = DMPOSum(rhoA, rhoB)
	% gather constants
	lengthA = size(rhoA, 1);
	lengthB = size(rhoB, 1);
	hilbyA = size(rhoA{1}, 3);
	hilbyB = size(rhoB{1}, 3);

	% check that rhoA and rhoB describe similar systems
	if lengthA ~= lengthB || hilbyA ~= hilbyB
		msgID = 'DMPOSum:BadArgs';
		msg = sprintf('Arguments must describe systems with the same local state dimensions, and the same number of sites. lengthA = %d, lengthB = %d, hilbyA = %d, hilbyB = %d', lengthA, lengthB, hilbyA, hilbyB);
		unmatchedDMPOSumArgs = MException(msgID, msg);		
		throw(unmatchedDMPOSumArgs);
	end

	HILBY = hilbyA;
	LENGTH = lengthA;

	% return allocation
	sumDMPO = cell(LENGTH, 1);

	% first and last sites are special snowflakes
	for ket = 1 : 1 : HILBY
		for bra = 1 : 1 : HILBY
			sumDMPO{1}(:, :, bra, ket) = horzcat(rhoA{1}(:, :, bra, ket), rhoB{1}(:, :, bra, ket));
			sumDMPO{LENGTH}(:, :, bra, ket) = vertcat(rhoA{LENGTH}(:, :, bra, ket), rhoB{LENGTH}(:, :, bra, ket));
		end
	end

	% and the rest...
	for site = 2 : 1 : LENGTH - 1
		for ket = 1 : 1 : HILBY
			for bra = 1 : 1 : HILBY
				sumDMPO{site}(:, :, bra, ket) = blkdiag(rhoA{site}(:, :, bra, ket), rhoB{site}(:, :, bra, ket));
			end
		end
	end
end
