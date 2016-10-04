% EigenSolver.m
% Accessor method for eigensolving routines which depend on whether the problem
% is using the Hermitian product or not, returns the largest real eigenvalue
% vector for the direct Liouvillian, and the smallest magnitude eigenvalue for
% the Hermitian product of the Liouvillian
% 2016-10-03

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN)
    if HERMITIAN
        opts.eps = 1E-12;
        [eigVector, eigValue] = primme_eigs(effL, 1, 'SA', opts);
    else
        opts.maxit = 500;
        opts.tol = eps;
        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
