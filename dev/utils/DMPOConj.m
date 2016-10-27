% DMPOConj.m
% function which returns the Hermitian conjugate of a density matrix
% product operator
% Oliver Thomson Brown
% 2016-02-08
%
% [ conjDMPO ] = DMPOConj(dmpo)
%
% RETURN
% conjDMPO:	cell array, the Hermitian conjugate of the input density
%			matrix product operator
%
% INPUT
% dmpo:		cell array, a density matrix product operator

function [conjDMPO] = DMPOConj(dmpo)
	% gather constants
	LENGTH = size(dmpo, 1);
	HILBY = size(dmpo{1}, 3);

	% allocate return
	conjDMPO = dmpo;

	for site = 1 : 1 : LENGTH
		for ket = 1 : 1 : HILBY
			for bra = 1 : 1 : HILBY
				conjDMPO{site}(:, :, ket, bra) = conj(dmpo{site}(:, :, bra, ket));
			end
		end
	end
end
