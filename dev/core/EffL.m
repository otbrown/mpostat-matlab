% EffL.m
% function which returns the effective Liouvillian for a particular site
% in the form of a matrix, ready to be eigensolved
% Oliver Thomson Brown
% 2016-05-05
%
% [effectiveLiouv] = EffL(TARGET, dmpo, mpo, left, right)
%
% RETURN
% effectiveLiouv    : complex double matrix, the effective Liovillian for the
%                     site specified by TARGET, reshaped to be a matrix
%
% INPUTS
% TARGET            : integer, the site for which the effective Liouvillian
%                     should be formed
% dmpo              : cell array, a density matrix product operator
% mpo               : cell array, mpo{n} contains the 6-dimensional matrix
%                     product operator tensor for the site n
% left              : cell array, contains the left blocks for each site
% right             : cell array, contains the right blocks for each site

function [effectiveLiouv] = EffL(TARGET, dmpo, mpo, left, right)
    % collate arguments for EffOpTensor
    [ROW_SIZE, COL_SIZE, HILBY, ~] = size(dmpo{TARGET});
    siteMPO = mpo{TARGET};
    [~, ~, ~, ~, OP_ROW, OP_COL] = size(siteMPO);
    lBlock = left{TARGET};
    rBlock = right{TARGET};

    % get effective operator from EffOpTensor
    effTensor = EffOpTensor(lBlock, siteMPO, rBlock, ROW_SIZE, COL_SIZE, ...
                            HILBY, OP_ROW, OP_COL);

    % reshape and return
    effTensor = permute(effTensor, [3, 4, 7, 8, 1, 2, 5, 6]);
    effectiveLiouv = reshape(effTensor, [ROW_SIZE*COL_SIZE*HILBY*HILBY, ...
                                         ROW_SIZE*COL_SIZE*HILBY*HILBY]);
end
