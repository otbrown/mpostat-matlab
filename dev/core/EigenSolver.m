% EigenSolver.m
% Interface function for eigensolving routines which depend on whether the
% problem is using the Hermitian product or not, returns the largest real
% eigenvalue vector for the direct Liouvillian, and the smallest magnitude
% eigenvalue vector for the Hermitian product of the Liouvillian
% Oliver Thomson Brown
% 2016-10-03
%
% [ eigVector, eigValue ] = EigenSolver(effL, HERMITIAN)
% [ eigVector, eigValue ] = EigenSolver(effL, HERMITIAN, initVec)
%
% RETURN
% eigVector:    (complex) double, the desired eigenvector of the effective %                Liouvillian
% eigValue:     (complex) double, the eigenvalue corresponding to eigVector
%
% INPUT
% effL:         (complex) double, effective Liouvillian for some site in
%               the system
% HERMITIAN:    bool, true if effL is hermitian, in which case PRIMME will
%               be used instead of eigs, as ARPACK struggles with smallest
%               magnitude eigenvalues
% initVec:      (complex) double, OPTIONAL, initial guess for the
%               eigenvector -- only used in the non-hermitian case PRIMME's
%               Matlab interface does not currently support this

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, varargin)
    if HERMITIAN
        % clean effL as primme has no tolerance for non-hermitian input
        epsilon = full(max(max(abs(effL - ctranspose(effL)))));
        fprintf('Hermiticity error: %g\n', epsilon);
        effL = (effL + ctranspose(effL))/2;

        opts = struct('eps', 1E-14);

        fprintf('\n');
        [eigVector, eigValue] = primme_eigs(effL, 1, 'SM', opts);
        fprintf('\n');
    else
        if nargin > 2
            opts.v0 = varargin{1};
        end

        opts.maxit = 500;
        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
