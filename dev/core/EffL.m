% EffL.m
% accessor method which returns the effective Liouvillian
% 2016-09-29
%
% [effectiveLiouv] = EffL(TARGET, dmpo, mpo, left, right)
%
% RETURN
% effectiveLiouv    : sparse complex double matrix, the effective Liovillian
%                     for the site specified by TARGET, reshaped to be a
%                     matrix
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
    % gather pass-forward variables
    [ROW_SIZE, COL_SIZE, HILBY, ~] = size(dmpo{TARGET});
    siteMPO = mpo{TARGET};
    lBlock = left{TARGET};
    rBlock = right{TARGET};

    effectiveLiouv = EffLSparse(lBlock, siteMPO, rBlock, ...
                                ROW_SIZE, COL_SIZE, HILBY);
end
