% DMPOScalarDiv.m
% divide a density matrix product operator by a scalar
% Oliver Thomson Brown
% 2016-02-08
%
% divDMPO = DMPOScalarDiv(dmpo, scalar)
%
% RETURN
% divDMPO:		cell array, rho / n
%
% INPUTS
% dmpo:			cell array, rho
% scalar:		double, n

function [divDMPO] = DMPOScalarDiv(dmpo, scalar)
	% will simply divide first tensor in the dmpo by
	% the scalar -- it is however possible to spread
	% the division out through the mpo
	% this may be desirable in order to prevent the first tensor
	% becoming too small, but is time costly and may lead
	% to more numerical errors depending on how far the
	% scalar division is spread

	divDMPO = dmpo;
	divDMPO{1} = divDMPO{1} / scalar;
end
