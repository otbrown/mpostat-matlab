% PhasedSearch.m
% top level function which carries out repeated stationary state searches,
% gradually lowering the accuracy threshold, and increasing the matrix
% dimension
% Oliver Thomson Brown
% 2017-02-07
%
% [ dmpoStat, phaseEigs ] = PhasedSearch(HILBY, LENGTH, mpo, ULTIMATE_THRESHOLD,
%                                         MAX_COMPRESS, VARIANT)
%
% RETURN
% dmpoStat:     cell, density matrix product operator representing the
%               stationary state (hopefully)
% phaseTrack:   (complex) double, contains the final eigenvalue of each
%               phase ofthe calculation
%
% INPUT
% HILBY:                integer, size of the local state space
% LENGTH:               integer, the number of sites in the system
% mpo:                  cell, Liouvillian for the system in matrix product
%                       operator form
% ULTIMATE_THRESHOLD:   double, how close must L*rho be to zero for the
%                       calculation to be deemed successful, the
%                       calculation will end once this is reached
% MAX_COMPRESS:         integer, the maximum allowed matrix dimension of
%                       the density matrix product operator
% VARIANT:              string, specifies whether to solve the
%                       non-Hermitian Liouvillian, or the Hermitian
%                       product, 'direct' or 'hermitian'

function [dmpoStat, phaseTrack] = PhasedSearch(HILBY, LENGTH, mpo, ULTIMATE_THRESHOLD, MAX_COMPRESS, VARIANT)
    % set up
    threshold = 0.01;
    compress = HILBY^2;
    phaseCount = 0;
    phaseTrack = [];
    phaseEig = Inf;
    eigTrack = NaN;
    ARPACK_msgID = 'MATLAB:eigs:ARPACKroutineErrorMinus14';

    dmpoStat = DDMPO(HILBY, LENGTH, compress);

    while phaseEig > ULTIMATE_THRESHOLD
        fprintf('PHASE %g:\n[ THRESHOLD: %g, COMPRESS: %g]\n', ...
                phaseCount, threshold, compress);

        dmpoStat = DMPOResize(dmpoStat, compress);

        try
            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(dmpoStat, mpo, ...
                                                threshold, VARIANT);
            phaseEig = abs(eigTrack(end));
        catch ME
            if strcmpi(ME.identifier, ARPACK_msgID)
                % DO NOTHING
            else
                throw(ME);
            end
        end

        % Technically this is dangerous as phaseTrack could grow to
        % infinite length. Best hope not, eh.
        phaseTrack(end + 1) = eigTrack(end);

        if phaseEig < threshold
            if phaseEig < ULTIMATE_THRESHOLD
                fprintf('PHASE SUCCESS. ULTIMATE THRESHOLD REACHED.\n');
                fprintf('ENDING CALCULATION\n');
            else
                fprintf('PHASE SUCCESS. LOWERING THRESHOLD.\n\n');
                threshold = max(ULTIMATE_THRESHOLD, threshold / 10);
            end
        else
            if compress == MAX_COMPRESS
                fprintf('PHASE FAILED. MATRIX DIMENSION LIMIT REACHED.\n');
                fprintf('ENDING CALCULATION.\n');
                break;
            else
                fprintf('PHASE FAILED. ENLARGING MATRIX DIMENSION.\n\n');
                compress = min(MAX_COMPRESS, compress + HILBY^2);
            end
        end

        phaseCount = phaseCount + 1;
    end
end
