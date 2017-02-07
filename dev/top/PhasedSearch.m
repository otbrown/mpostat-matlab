% PhasedSearch.m
% top level function which carries out repeated stationary state searches,
% gradually lowering the accuracy threshold, and increasing the matrix dimension
% Oliver Thomson Brown
% 2017-02-07
%
% [ dmpoStat, phaseEigs ] = PhasedSearch(dmpoInit, mpo, ULTIMATE_THRESHOLD,
%                                         MAX_COMPRESS, variant)
%
% RETURN
% dmpoStat:   cell, density matrix product operator representing the stationary
%             state (hopefully)
% phaseEigs:  (complex) double, contains the final eigenvalue of each phase of
%             the calculation
%
% INPUT
% HILBY:              integer, size of the local state space
% LENGTH:             integer, the number of sites in the system
% mpo:                cell, Liouvillian for the system in matrix product
%                     operator form
% ULTIMATE_THRESHOLD: double, how close must L*rho be to zero for the
%                     calculation to be deemed successful, the calculation will
%                     end once this is reached
% MAX_COMPRESS:       integer, the maximum allowed matrix dimension of the
%                     density matrix product operator
% variant:            string, specifies whether to solve the non-Hermitian
%                     Liouvillian, or the Hermitian product, 'direct' or
%                     'hermitian'
