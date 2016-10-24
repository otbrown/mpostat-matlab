% EigenSolver.m
% Accessor method for eigensolving routines which depend on whether the problem
% is using the Hermitian product or not, returns the largest real eigenvalue
% vector for the direct Liouvillian, and the smallest magnitude eigenvalue for
% the Hermitian product of the Liouvillian
% 2016-10-03

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN, varargin)
    if HERMITIAN
        % clean effL as primme has no tolerance for non-hermitian input
        epsilon = full(max(max(abs(effL - ctranspose(effL)))));
        fprintf('Hermiticity error: %g\n', epsilon);
        effL = (effL + ctranspose(effL))/2;

        opts = struct('numTargetShifts', 1, 'targetShifts', [0], 'eps', 1E-14);

        fprintf('\n');
        [eigVector, eigValue] = primme_eigs(effL, 1, 'CT', opts);
        fprintf('\n');
    else
        if nargin > 2
            opts.v0 = varargin{1};
        end

        opts.maxit = 500;
        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
