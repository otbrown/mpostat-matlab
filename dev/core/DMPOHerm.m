% DMPOHerm.m
% makes the supplied density matrix product operator Hermitian, by
% calculating (rho + ctranspose(rho)) / 2
% Oliver Thomson Brown
% 2016-02-17
%
% [ hermDMPO ] = DMPOHerm(dmpo)
%
% RETURN
% hermDMPO: cell, a density matrix product operator corresponding to
%           (rho + ctranspose(rho)) / 2
%
% INPUT
% dmpo: cell, a density matrix product operator representing the density
%       matrix rho

function hermDMPO = DMPOHerm(dmpo)
    hermDMPO = DMPOScalarDiv(DMPOSum(dmpo, DMPOConj(dmpo)), 2);
end
