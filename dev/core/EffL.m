% EffL.m
% accessor method which returns the effective Liouvillian built by either
% EffLFull or EffLSparse, depending on the value of logical MEMSAVE
% 2016-09-29
%
% [effectiveLiouv] = EffL(TARGET, dmpo, mpo, left, right, MEMSAVE)
%
% RETURN
% effectiveLiouv    : complex double matrix, the effective Liovillian for the
%                     site specified by TARGET, reshaped to be a matrix, may be
%                     full or space, depending on MEMSAVE
%
% INPUTS
% TARGET            : integer, the site for which the effective Liouvillian
%                     should be formed
% dmpo              : cell array, a density matrix product operator
% mpo               : cell array, mpo{n} contains the 6-dimensional matrix
%                     product operator tensor for the site n
% left              : cell array, contains the left blocks for each site
% right             : cell array, contains the right blocks for each site
% MEMSAVE           : if true EffL passes on input to EffLSparse, otherwise it
%                     calls EffLFull

function [effectiveLiouv] = EffL(TARGET, dmpo, mpo, left, right, MEMSAVE)
    % gather pass-forward variables
    [ROW_SIZE, COL_SIZE, HILBY, ~] = size(dmpo{TARGET});
    siteMPO = mpo{TARGET};
    lBlock = left{TARGET};
    rBlock = right{TARGET};

    if MEMSAVE
        effectiveLiouv = EffLSparse(lBlock, siteMPO, rBlock, ...
                                    ROW_SIZE, COL_SIZE, HILBY);
    else
        effectiveLiouv = EffLFull(lBlock, siteMPO, rBlock, ...
                                  ROW_SIZE, COL_SIZE, HILBY);
    end
end
