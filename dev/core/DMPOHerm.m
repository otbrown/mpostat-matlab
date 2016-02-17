% DMPOHerm.m
% makes the supplied density matrix product operator Hermitian, by calculating
% (rho + ctranspose(rho)) / 2
% Oliver Thomson Brown
% 2016-02-17
%
% hermDMPO = DMPOHerm(dmpo)
%
% RETURN
% hermDMPO      : a density matrix product operator corresponding to
%                 (dmpo + ctranspose(dmpo)) / 2
%
% INPUTS
% dmpo          : an arbitrary density matrix product operator

function hermDMPO = DMPOHerm(dmpo)
    hermDMPO = DMPOScalarDiv(DMPOSum(dmpo, DMPOConj(dmpo)), 2);
end
