% TrNorm.m
% trace normalises an arbitrary density matrix product operator
% Oliver Thomson Brown
% 2016-02-17
%
% normDMPO = TrNorm(dmpo)
%
% RETURN
% normDMPO  : cell array, contains the trace normalised dmpo
%
% INPUTS
% dmpo      : cell array, an arbitrary dmpo

function [normDMPO] = TrNorm(dmpo)
    tr = DMPOTrace(dmpo);
    normDMPO = DMPOScalarDiv(dmpo, tr);
end
