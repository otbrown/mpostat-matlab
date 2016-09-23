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

function [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD)
    % gather up system parameters
    LENGTH = length(dmpoInit);
    HILBY = size(dmpoInit{1}, 3);
    MIDDLE = floor(LENGTH/2);
    MAX_DIM = size(dmpoInit{MIDDLE}, 2);
    MAX_LDIM = (MAX_DIM * HILBY)^4;

    % set internal calculation parameters
    RUNMAX = 50*LENGTH;
    CONVERGENCE_THRESHOLD = THRESHOLD / (2 * LENGTH);
    convFlag = false;
    success = false;
    finished = false;
    sweepCount = 0;
    updCount = 0;
    opts.maxit = 500;
    route = 1 : 1 : LENGTH;

    % print some info about the calculation
    fprintf('Variational Stationary State Search\n');
    fprintf('%s\n\n', datestr(datetime('now'), 31));
    fprintf('System Parameters:\n\tSystem size: %g\n\tLocal states: %g\n', LENGTH, HILBY);
    fprintf('Calculation Parameters:\n\tEigenvalue threshold: %g\n\tMaximum MPS matrix size: %g\n\tMaximum effective Liouvillian size: %g\n\n', THRESHOLD, MAX_DIM, MAX_LDIM);

    % allocate return variables
    % CLEAN UP INITIAL STATE
    dmpoStat = dmpoInit;
    eigTrack = NaN(LENGTH, 1);

    % build left and right blocks
    % MAKE AN INTERFACE FUNCTION LIKE GROWBLOCK
    left = cell(LENGTH,1);
    right = cell(LENGTH, 1);
    left{1} = 1;
    right{LENGTH} = 1;

    for site = LENGTH : -1 : 2
        [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
        [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
        dmpoStat = RCan(dmpoStat, site);
        right{site - 1} = GrowRight(dmpoStat{site}, mpo{site}, right{site}, ...
                                    ROW_SIZE, COL_SIZE, HILBY, OP_ROW);
    end

    while ~convFlag && updCount < RUNMAX
        for site = route
            effL = EffL(site, dmpoStat, mpo, left, right);
            [update, eigTrack(LENGTH)] = eigs(effL, 1, 'lr', opts);

            [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});

            for bra = 0 : 1 : (HILBY - 1)
                for ket = 0 : 1 : (HILBY - 1)
                    for row = 0 : 1 : (ROW_SIZE - 1)
                        for col = 1 : 1 : COL_SIZE
                            dmpoStat{site}(row+1, col, bra+1, ket+1) = update(ket*HILBY*ROW_SIZE*COL_SIZE + bra*ROW_SIZE*COL_SIZE + row*COL_SIZE + col);
                        end
                    end
                end
            end

            % canonicalise & include new site in block tensor
            if mod(sweepCount, 2)
                if site ~= 1
                    % RCAN
                    dmpoStat = RCan(dmpoStat, site);
                    [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
                    [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
                    right{site - 1} = GrowRight(dmpoStat{site}, mpo{site}, right{site}, ...
                                         ROW_SIZE, COL_SIZE, HILBY, OP_ROW);
                end
            else
                if site ~= LENGTH
                % LCAN
                dmpoStat = LCan(dmpoStat, site);
                [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
                [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpo{site});
                left{site + 1} = GrowLeft(dmpoStat{site}, mpo{site}, left{site}, ...
                                         ROW_SIZE, COL_SIZE, HILBY, OP_COL);
                end
            end

            % evaluate convergence
            [convFlag, convergence] = ConvTest(eigTrack, CONVERGENCE_THRESHOLD);

            % stop following route if RUNMAX is reached or if the calculation
            % has converged to desired threshold
            if convFlag || updCount == RUNMAX
                finished = true;
                if abs(eigTrack(LENGTH)) < THRESHOLD
                    fprintf('\nCalculation successful.\n');
                    success = true;
                    break;
                else
                    fprintf('\nCalculation failed to reach desired accuracy. Larger matrix dimensions may be required.\n');
                    break;
                end
            end

            % drop oldest eigenvalue from eigTrack and move elements back
            eigTrack(1 : LENGTH - 1) = eigTrack(2 : LENGTH);
            updCount = updCount + 1;
        end

        if finished
            fprintf('Finished at %s.\n[ Eigenvalue: %g, Convergence: %g ]\n', datestr(datetime('now'), 31), eigTrack(LENGTH), convergence);
        else
            % add to sweepCount and report on progress
            sweepCount = sweepCount + 1;
            fprintf('Sweep %g:\n[ Eigenvalue: %g, Convergence: %g ]\n', sweepCount, eigTrack(LENGTH), convergence);
        end

        % flip it and reverse it
        route = flip(route);
    end
    % trace normalise the state
    dmpoStat = TrNorm(dmpoStat);
end
