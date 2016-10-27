% TrNorm.m
% trace normalises an arbitrary density matrix product operator
% Oliver Thomson Brown
% 2016-02-17
%
% [ normDMPO ] = TrNorm(dmpo)
%
% RETURN
% normDMPO: cell, contains the trace normalised dmpo
%
% INPUT
% dmpo: cell, a density matrix product operator

function [normDMPO] = TrNorm(dmpo)
    tr = DMPOTrace(dmpo);
    normDMPO = DMPOScalarDiv(dmpo, tr);
end
