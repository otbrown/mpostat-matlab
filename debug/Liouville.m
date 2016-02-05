% Liouville.m
% function which generates the Liouvillian matrix from the Hamiltonian of a
% system, and the dissipative terms
% Oliver Thomson Brown
% 2015-06-10
%
% Return
% liouvillian   : 2D array, containing the coefficient matrix L for the
%                   equation L*rho = drho/dt
%
% Inputs
% hamiltonian   : hamiltonian for the system
% dissipators   : cell array, 1 x numDiss x 2, dissipators{1, dissdex, 1}
%                   contains a dissipative constant, and dissipators{1,
%                   dissdex, 2} contains the operator O, where the first 
%                   term in the dissipative part of the equation is 
%                   2*O*rho*O^{\dagger} 
% HILBY         : size of the local Hilbert space
% NUM_SITES     : number of sites in the system

function [liouvillian] = Liouville(hamiltonian, dissipators, HILBY, NUM_SITES)
    % constants
    dimension = HILBY^NUM_SITES;
    dissNum = size(dissipators, 2);
    
    % pralloc
    densityMatrix = sym('rho', [dimension, dimension]);
    rhoVec = reshape(densityMatrix, [1, dimension^2]);
    lindDiss = zeros(dimension);
    
    % commutator
    lindComm = -1i * (hamiltonian * densityMatrix - densityMatrix * hamiltonian);
    
    % dissipative terms
    for dissdex = 1 : 1 : dissNum
        dissConstant = dissipators{1, dissdex, 1};
        dissOp = dissipators{1, dissdex, 2};
        conjDissOp = ctranspose(dissOp);
        lindDiss = lindDiss + dissConstant * ...
            ( 2 * dissOp * densityMatrix * conjDissOp ...
            - conjDissOp * dissOp * densityMatrix ...
            - densityMatrix * conjDissOp * dissOp ) / 2;
    end
    
    % form drho/dt and reshape
    lindblad = lindComm + lindDiss;
    dRhoVec = reshape(lindblad, [dimension^2, 1]);
    
    % generate liouvillian
    liouvillian = equationsToMatrix(dRhoVec, rhoVec);
end