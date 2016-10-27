% DMPOScalarDiv.m
% divide a density matrix product operator by a scalar
% Oliver Thomson Brown
% 2016-02-08
%
% [ divDMPO ] = DMPOScalarDiv(dmpo, scalar)
%
% RETURN
% divDMPO:	cell array, a dmpo -- rho / n
%
% INPUT
% dmpo:		cell array, a density matrix product operator -- rho
% scalar:	(complex) double, n

function [divDMPO] = DMPOScalarDiv(dmpo, scalar)
	divDMPO = dmpo;
	divDMPO{1} = divDMPO{1} / scalar;
end
