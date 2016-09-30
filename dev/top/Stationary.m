% Stationary.m
% top level function for the stationary state search
% Oliver Thomson Brown
% 2016-05-06
%
% [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD)
% [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, variant)
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
% variant       : string, optional argument, may specify whether to build the
%                 effective Liouvillian in full or as a sparse array -- valid
%                 strings are 'full' and 'sparse' -- 'sparse' is default

function [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, varargin)
    % check for optional arguments
    switch nargin
        case 3
            MEMSAVE = true;
        case 4
            if strcmpi(varargin{1}, 'full')
                MEMSAVE = false;
            elseif strcmpi(varargin{1}, 'sparse')
                MEMSAVE = true;
            else
                ME = MException('Stationary:badMEMSAVE', 'The last argument was invalid: %s. Type help Stationary.', varargin{1});
                throw(ME);
            end
        otherwise
            ME = MException('Stationary:badArguments', 'Stationary accepts 3 or 4 arguments, but %g were supplied. Type help Stationary.', nargin);
            throw(ME);
    end

    % gather up physical system parameters
    LENGTH = length(dmpoInit);
    HILBY = size(dmpoInit{1}, 3);
    MIDDLE = floor(LENGTH/2);
    MAX_DIM = size(dmpoInit{MIDDLE}, 2);
    MAX_LDIM = (MAX_DIM * HILBY)^2;

    % set internal calculation parameters
    RUNMAX = 50*LENGTH;
    CONVERGENCE_THRESHOLD = THRESHOLD / (2 * LENGTH);

    % set some parameters for 'finessing' eigs
    opts.maxit = 500;

    % initialise flags and counters
    convFlag = false;
    success = false;
    finished = false;
    sweepCount = 0;
    updCount = 0;
    route = 1 : 1 : LENGTH;

    % make mpo hermitian
    mpo = MPOHermProd(mpo);

    % print some info about the calculation
    fprintf('Variational Stationary State Search\n');
    fprintf('%s\n\n', datestr(datetime('now'), 31));
    fprintf('System Parameters:\n\tSystem size: %g\n\tLocal states: %g\n', LENGTH, HILBY);
    fprintf('Calculation Parameters:\n\tEigenvalue threshold: %g\n\tMaximum MPS matrix size: %g\n\tMaximum effective Liouvillian size: %g\n', THRESHOLD, MAX_DIM, MAX_LDIM);
    if MEMSAVE
        fprintf('\tEffective Liouvillian representation: Sparse\n\n');
    else
        fprintf('\tEffective Liouvillian representation: Full\n\n');
    end

    % allocate return variables
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
            effL = EffL(site, dmpoStat, mpo, left, right, MEMSAVE);
            [update, eigTrack(LENGTH)] = eigs(effL, 1, 'sm', opts);

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
