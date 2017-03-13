% Stationary.m
% top level function for the stationary state search
% Oliver Thomson Brown
% 2016-05-06
%
% [ dmpoStat, eigTrack ] = Stationary(dmpoInit, mpo, THRESHOLD)
% [ dmpoStat, eigTrack ] = Stationary(dmpoInit, mpo, THRESHOLD, variant)
%
% RETURN
% dmpoStat: cell, density matrix product operator representing the
%           stationary state (hopefully)
% eigTrack: (complex) double, contains the eigenvalues from the last (2 *
%           (LENGTH-1)) site updates which are monitored for convergence
%
% INPUT
% dmpoInit:     cell, a density matrix product operator representing some
%               initial state
% mpo:          cell, Liouvillian for the system in matrix product
%               operator form
% THRESHOLD:    double, how close must L*rho be to zero for the
%               calculation to be deemed successful
% variant:      string, specify whether to solve the non-Hermitian
%               Liovillian, or the Hermitian product L^(T*)L --
%               'hermitian' and 'direct' are the two accepted values

function [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, variant)

    if strcmpi(variant, 'hermitian')
        HERMITIAN = true;
    elseif strcmpi(variant, 'direct')
        HERMITIAN = false;
    else
        ME = MException('Stationary:badHERMITIAN', ['The last ',...
        'argument was invalid: %s. Type help Stationary.'], ...
        variant);
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

    % print some info about the calculation
    fprintf('Variational Stationary State Search\n');
    fprintf('%s\n\n', datestr(datetime('now'), 31));
    fprintf(['System Parameters:\n\tSystem size: %g\n\tLocal ', ...
            'states: %g\n'], LENGTH, HILBY);
    fprintf(['Calculation Parameters:\n\tEigenvalue threshold: ', ...
            '%g\n\tMaximum MPS matrix size: %g\n\tMaximum effective ', ...
            'Liouvillian size: %g\n'], THRESHOLD, MAX_DIM, MAX_LDIM);
    if HERMITIAN
        fprintf('\tEffective Liouvillian: Hermitian Product\n\n');
    else
        fprintf('\tEffective Liouvillian: Non-Hermitian\n\n');
        ARPACK_msgID = 'MATLAB:eigs:ARPACKroutineErrorMinus14';
    end

    % allocate return variables
    dmpoStat = dmpoInit;
    eigTrack = NaN(2 * (LENGTH-1), 1);

    % build left and right blocks
    left = cell(LENGTH,1);
    right = cell(LENGTH, 1);
    left{1} = 1;
    right{LENGTH} = 1;

    for site = LENGTH : -1 : 2
        direction = 'R';
        dmpoStat = Can(dmpoStat, site, direction);
        right{site - 1} = GrowBlock(dmpoStat, mpo, left, right, site, ...
                                    direction);
    end

    % initialise flags and counters
    convFlag = false;
    finished = false;
    sweepCount = 0;
    updCount = 0;
    route = 1 : 1 : LENGTH;
    direction = 'L';

    while ~convFlag && updCount < RUNMAX
        for sitedex = 1 : 1 : (LENGTH - 1)
            site = route(sitedex);

            [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpoStat{site});
            effL = EffL(site, dmpoStat, mpo, left, right);

            if HERMITIAN
                [update, eig] = EigenSolver(effL, HERMITIAN, THRESHOLD);
            else
                % we can supply an initial guess to eigs, so we use the
                % current site tensor, to aid convergence
                siteVec = permute(dmpoStat{site}, [2, 1, 3, 4]);
                siteVec = reshape(siteVec, [ROW_SIZE*COL_SIZE*HILBY^2, 1]);
                try
                    [update, eig] = EigenSolver(effL, HERMITIAN, ...
                                                THRESHOLD, siteVec);
                catch ME
                    if strcmp(ME.identifier, ARPACK_msgID)
                        fname = sprintf('mpostat%uX%u_failed.mat', ...
                                        LENGTH, HILBY);
                        save(fname);
                        fprintf(['Unfortunately, the calculation has ' ...
                                'failed while trying to find  ' ...
                                'eigenvalues.\nPartial result saved in %s\nConsider using larger matrix ' ...
                                'dimensions or the Hermitian ' ...
                                'variant.\n'], fname);
                        throw(ME);
                    else
                        throw(ME);
                    end
                end
            end

            if ~(abs(eig) > abs(eigTrack(end)))
                eigTrack(end) = eig;
                update = reshape(update, ...
                                    [COL_SIZE, ROW_SIZE, HILBY, HILBY]);
                dmpoStat{site} = permute(update, [2, 1, 3, 4]);
            end

            % canonicalise & include new site in block tensor
            dmpoStat = Can(dmpoStat, site, direction);

            if direction == 'L'
                left{site + 1} = GrowBlock(dmpoStat, mpo, left, right, ...
                                            site, direction);
            else
                right{site - 1} = GrowBlock(dmpoStat, mpo, left, right, ...
                                            site, direction);
            end

            % evaluate convergence
            [convFlag, convergence] = ...
                ConvTest(eigTrack, CONVERGENCE_THRESHOLD);

            % stop following route if RUNMAX is reached or if the
            % calculation has converged to desired threshold
            if convFlag || updCount == RUNMAX
                finished = true;
                if abs(eigTrack(end)) < THRESHOLD
                    fprintf('\nCalculation successful.\n');
                    break;
                else
                    fprintf(['\nCalculation failed to reach desired ', ...
                            'accuracy. Larger matrix dimensions may ', ...
                            'be required.\n']);
                    break;
                end
            end

            % drop oldest eigenvalue from eigTrack and move elements back
            eigTrack(1 : (end - 1)) = eigTrack(2 : end);
            updCount = updCount + 1;
        end

        if finished
            fprintf(['Finished at %s.\n', ...
                    '[ Eigenvalue: %g, Convergence: %g ]\n'], ...
                    datestr(datetime('now'), 31), eigTrack(end), ...
                    convergence);
        else
            % add to sweepCount and report on progress
            sweepCount = sweepCount + 1;
            fprintf('Sweep %g:\n[ Eigenvalue: %g, Convergence: %g ]\n', ...
                    sweepCount, eigTrack(end), convergence);
        end

        % flip it and reverse it
        route = flip(route);
        if direction == 'L'
            direction = 'R';
        else
            direction = 'L';
        end
    end
    % trace normalise the state
    dmpoStat = TrNorm(dmpoStat);
end
