% Stationary.m
% top level function for the stationary state search
% Oliver Thomson Brown
% 2016-05-06
%
% [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, RUNMAX)
%
% RETURN
% dmpoStat      : cell array, density matrix product operator representing the
%                 stationary state (hopefully)
% eigTrack      : double array, contains the eigenvalues from each site update
%                 monitored for convergence
%
% INPUTS
% dmpoInit      : cell array, density matrix product operator containing some
%                 initial state
% mpo           : cell array, Liouvillian for the system in matrix product
%                 operator form
% THRESHOLD     : double, the convergence threshold
% RUNMAX        : integer, the maximum number of updates before the code fails

function [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, RUNMAX)
    % gather up system parameters
    LENGTH = length(dmpoInit);
    HILBY = size(dmpoInit{1}, 3);

    % allocate return variables
    % CLEAN UP INITIAL STATE
    dmpoStat = dmpoInit;
    eigTrack = NaN(RUNMAX - 1, 1);

    % build left and right blocks
    % MAKE AN INTERFACE FUNCTION LIKE GROWBLOCK
    left = cell(LENGTH,1);
    right = cell(LENGTH, 1);
    left{1} = 1;
    right{LENGTH} = 1;

    for site = 1 : 1 : (LENGTH - 1)
        [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
        [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
        left{site + 1} = GrowLeft(dmpoStat{site}, mpo{site}, left{site}, ...
                                  ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL);
    end

    for site = LENGTH : -1 : 2
        [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
        [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
        right{site - 1} = GrowRight(dmpoStat{site}, mpo{site}, right{site}, ...
                                    ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL);
    end

    % Run the search
    convFlag = 0;
    sweepCount = 0;
    updCount = 1;
    route = 1 : 1 : LENGTH;

    while ~convFlag && updCount < RUNMAX
        for site = route
            [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
            [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
            effectiveLiouv = EffL(site, dmpoStat, mpo, left, right);
            effectiveLiouv = effectiveLiouv + ...
                             eps*eye(HILBY*HILBY*ROW_SIZE*COL_SIZE);
            opts.tol = eps;
            [update , eigTrack(updCount)] = eigs(effectiveLiouv, 1, 'sm', opts);

            dmpoStat{site} = reshape(update,[ROW_SIZE, COL_SIZE, HILBY, HILBY]);

            updCount = updCount + 1;

            % canonicalise & include new site in block tensor
            if mod(sweepCount, 2) && site ~= 1;
                % RCAN
                dmpoStat = RCan(dmpoStat, site);
                [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
                right{site - 1} = GrowRight(dmpoStat{site}, mpo{site}, ...
                                            right{site}, ROW_SIZE, COL_SIZE, ...
                                            HILBY, OP_ROW, OP_COL);
            elseif site ~= LENGTH
                % LCAN
                dmpoStat = LCan(dmpoStat, site);
                [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
                left{site + 1} = GrowLeft(dmpoStat{site}, mpo{site}, ...
                                          left{site}, ROW_SIZE, COL_SIZE, ...
                                          HILBY, OP_ROW, OP_COL);
            end

            % ADD CONVERGENCE TEST!

            % stop following route if RUNMAX is reached
            if updCount == RUNMAX
                dmpoStat = TrNorm(dmpoStat);
                break;
            end
        end

        dmpoStat = TrNorm(dmpoStat);
        sweepCount = sweepCount + 1;
        % flip it and reverse it
        route = flip(route);
    end
end
