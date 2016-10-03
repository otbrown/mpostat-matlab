% EigenSolver.m
% Accessor method for eigensolving routines which depend on whether the problem
% is using the Hermitian product or not, returns the largest real eigenvalue
% vector for the direct Liouvillian, and the smallest magnitude eigenvalue for
% the Hermitian product of the Liouvillian
% 2016-10-03

function [eigVector, eigValue] = EigenSolver(effL, HERMITIAN)
    opts.maxit = 500;

    if HERMITIAN
        [eigVector, eigValue] = eigs(effL, 1, 'sm', opts);
    else
        [eigVector, eigValue] = eigs(effL, 1, 'lr', opts);
    end
end
