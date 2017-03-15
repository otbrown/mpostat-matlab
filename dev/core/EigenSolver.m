% EigenSolver.m
% Interface function for eigensolving routines which depend on whether the
% problem is using the Hermitian product or not, returns the largest real
% eigenvalue vector for the direct Liouvillian, and the smallest magnitude
% eigenvalue vector for the Hermitian product of the Liouvillian
% Oliver Thomson Brown
% 2016-10-03
%
% [eigVector, eigValue] = EigenSolver(effL, HERMITIAN = true,
%                                       HERMITICITY_THRESHOLD)
% [eigVector, eigValue] = EigenSolver(effL, HERMITIAN = false, initVec)
%
% RETURN
% eigVector:    (complex) double, the desired eigenvector of the effective
%                Liouvillian
% eigValue:     (complex) double, the eigenvalue corresponding to eigVector
%
% INPUT
% effL:                     (complex) double, effective Liouvillian for
%                           some site in the system
% HERMITIAN:                bool, true if effL is hermitian, in which case
%                           PRIMME will be used instead of eigs, as ARPACK
%                           struggles with smallest magnitude eigenvalues
% HERMITICITY_THRESHOLD:    (complex) double, OPTIONAL, should ONLY be
%                           supplied if HERMITIAN = true, threshold at
%                           which function will reject mpo as probably not
%                           actually Hermitian -- in Stationary this is
%                           set at 0.1 * the smallest value in the mpo
% initVec:                  (complex) double, OPTIONAL, initial guess for
%                           the eigenvector -- only used in the
%                           non-hermitian case PRIMME's Matlab interface
%                           does not currently support this

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, initVec, varargin)
    narginchk(3, 4);

    if HERMITIAN
        % solve hermitian product of L
        if nargin == 4
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

        opts = struct('eps', 1E-14, 'v0', {initVec});

        [eigVector, eigValue] = primme_eigs(effL, 1, 'SM', opts);
    else
        % solve L directly
        opts = struct('maxit', 500, 'v0', initVec);
        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
