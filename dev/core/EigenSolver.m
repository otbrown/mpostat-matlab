% EigenSolver.m
% Interface function for eigensolving routines which depend on whether the
% problem is using the Hermitian product or not, returns the largest real
% eigenvalue vector for the direct Liouvillian, and the smallest magnitude
% eigenvalue vector for the Hermitian product of the Liouvillian
% Oliver Thomson Brown
% 2016-10-03
%
% [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, PRIMME, initVec)
% [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, PRIMME, initVec,
%                                       HERMITICITY_THRESHOLD)
%
% RETURN
% eigVector:    (complex) double, the desired eigenvector of the effective
%                Liouvillian
% eigValue:     (complex) double, the eigenvalue corresponding to eigVector
%
% INPUT
% effL:                     (complex) double, effective Liouvillian for
%                           some site in the system
% HERMITIAN:                bool, true if effL is hermitian
% PRIMME:                   bool, true if PRIMME eigensolver should be used
%                           instead of eigs, which can struggle with
%                           finding the smallest magnitude eigenvalue
% initVec:                  (complex) double, initial guess for the
%                           eigenvector -- only used with eigs PRIMME's
%                           Matlab interface does not currently support
%                           this
% HERMITICITY_THRESHOLD:    (complex) double, OPTIONAL, should ONLY be
%                           supplied if HERMITIAN = true, threshold at
%                           which function will reject mpo as probably not
%                           actually Hermitian -- in Stationary this is
%                           set at 0.1 * the smallest value in the mpo

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, PRIMME, initVec, varargin)
    narginchk(4, 5);

    if HERMITIAN
        % solve hermitian product of L
        if nargin == 5
            % if HERMITICITY_THRESHOLD is supplied check difference
            % between L and L'
            epsilon = full(max(max(abs(effL - ctranspose(effL)))));

            if epsilon > varargin{1}
                ME = MException('EigenSolver:badHermiticity',  ...
                ['The error in L'' - L was large. Supplied MPO', ...
                 ' may not be Hermitian.']);
                throw(ME);
            end
        end

        % clean effL as primme has no tolerance for non-hermitian input
        effL = (effL + ctranspose(effL))/2;

        if PRIMME
            opts = struct('eps', 1E-14, 'numTargetShifts', {1}, ...
                    'targetShifts', {0}, 'initialevecs', {initVec});

            fprintf('\n');
            [eigVector, eigValue] = primme_eigs(effL, 1, 'CT', opts);
            fprintf('\n');
        else
            opts = struct('maxit', 500, 'v0', initVec);

            [eigVector, eigValue] = eigs(effL, 1, 'sm', opts);
        end
    else
        % solve L directly
        opts = struct('maxit', 500, 'v0', initVec);

        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
