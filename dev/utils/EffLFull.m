% EffLFull.m
% function which returns the effective Liouvillian for a particular site
% in the form of a matrix, ready to be eigensolved
% Oliver Thomson Brown
% 2016-05-05
%
% [effectiveLiouv] = EffLFull(TARGET, dmpo, mpo, left, right)
%
% RETURN
% effectiveLiouv    : complex double matrix, the effective Liovillian for the
%                     site specified by TARGET, reshaped to be a matrix
%
% INPUTS
% lBlock          : 3-dimensional complex double array, rank-3 left block tensor
%                   which contains the contraction through the system from the
%                   first site up to the target site
% siteMPO         : 6-dimensional complex double array, rank-6 matrix product
%                   operator tensor for the target site
% rBlock          : 3-dimensional complex double array, rank-3 right block
%                   right block tensor which contains the contraction through
%                   the sytem from the target site to the last
% ROW_SIZE        : double, the size of the first virtual dimension of the
%                   density matrix product operator for the target site
% COL_SIZE        : double, the size of the second virtual dimension of the
%                   density matrix product operator for the target site
% HILBY           : double, the size of the local state space

function [effectiveLiouv] = EffLFull(lBlock, siteMPO, rBlock, ROW_SIZE, COL_SIZE, HILBY)
    % get effective operator from EffOpTensor
    effTensor = EffOpTensor(lBlock, siteMPO, rBlock, ROW_SIZE, COL_SIZE, HILBY);

    % reshape and return
    effTensor = permute(effTensor, [3, 4, 7, 8, 2, 1, 5, 6]);
    effectiveLiouv = reshape(effTensor, [ROW_SIZE*COL_SIZE*HILBY*HILBY, ...
                                         ROW_SIZE*COL_SIZE*HILBY*HILBY]);
end
