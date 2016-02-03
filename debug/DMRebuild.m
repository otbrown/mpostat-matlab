% DMRebuild.m
% a function to rebuild a density matrix from a density matrix MPS
% Oliver Thomson Brown
% 2015-08-07
%
% Valid calls:
% densityMatrix = DMRebuild(dmps)
% densityMatrix = DMRebuild(dmps, matrixFlag)
% [densityMatrix, epsilon] = DMRebuild(dmps, matrixFlag, actual) 
%
% Returns:
% !REQUIRED
% densityMatrix : dim(H)^2N vector, or dim(H)^N square array, complex double -- contains the density matrix recreated from the MPS
% !OPTIONAL
% epsilon       : dim(H)^2N vector, or dim(H)^N square array, complex double -- contains the absolute difference between the density matrix recreated from the MPS and another supplied
%
% Inputs:
% !REQUIRED
% dmps          : cell array, N by 1, where N is the length of the system, contains the density matrix mps matrices
% !OPTIONAL
% matrixFlag    : bool, jk double bc MATLAB. if ~= 0 then densityMatrix is returned as a square matrix, rather than a column vector
% actual        : dim(H)^2N vector, or dim(H)^N square array, complex double -- the density matrix against which the rebuilt one should be compared

function [densityMatrix, varargout] = DMRebuild(dmps, varargin)
    % INPUT AND RETURN CHECKS AND LABELLING
    nargoutchk(0,2);
    narginchk(1,3);

    LENGTH = size(dmps,1);
    HILBY = size(dmps{1}, 3);
    SPACE = HILBY^LENGTH;

    matrixFlag = 0;

    if nargout == 2
        if nargin < 3
            fprintf('You must provide the expected density matrix or vector, if you wish a difference matrix or vector to be returned.\nPlease also note that this should be supplied as the third input, the second input should be a 0 if you want to build the vectorised density matrix, or a 1 if you want it returned as a matrix.\n');
            return;
        else
            matrixFlag = varargin{1};
            actual = varargin{2};
            assert(isscalar(matrixFlag), 'matrixFlag (the second input argument) should be a 0 if you want the vectorised density matrix to be returned, or a 1 if you want it in matrix form');
            assert(numel(actual) == SPACE^2, 'actual (the third input argument) appears to be a different size from the density matrix which will be rebuilt');
            epsilon = zeros(SPACE^2, 1);
        end
        diffFlag = 1;
    else
        if nargin > 2
            fprintf('You must provide a return variable for the difference array: [densityMatrix, epsilon] = DMRebuild(dmps, matrixFlag, actual)\n');
            return;
        elseif nargin == 2
            matrixFlag = varargin{1};
            assert(isscalar(matrixFlag), 'matrixFlag (the second input argument) should be a 0 if you want the vectorised density matrix to be returned, or a 1 if you want it in matrix form');
        end
        diffFlag = 0;
    end

    % RETURN ALLOCATION
    densityMatrix = complex(zeros(SPACE^2, 1));
    
    % REBUILD VECTORISED DENSITY MATRIX
    for ketState = 0 : 1 : (SPACE - 1)
        for braState = 0 : 1 : (SPACE - 1)
            stateDex = braState * SPACE + ketState + 1;
            braBits = FWBase(braState, HILBY, LENGTH);
            ketBits = FWBase(ketState, HILBY, LENGTH);

            coefft = 1;
            for site = 1 : 1 : LENGTH
                bra = braBits(site) + 1;
                ket = ketBits(site) + 1;
                coefft = coefft * dmps{site}(:, :, bra, ket);
            end

            densityMatrix(stateDex) = coefft;

        end
    end

    % PREPARE OUTPUTS
    % difference array
    if diffFlag
        if size(actual, 2) ~= 1
            actual = reshape(actual, [SPACE^2, 1]);
        end
        epsilon = abs(actual) - abs(densityMatrix);
        if matrixFlag
            varargout{1} = reshape(epsilon, [SPACE, SPACE]);
        else
            varargout{1} = epsilon;
        end
    end
    % density matrix
    if matrixFlag
        densityMatrix = reshape(densityMatrix, [SPACE, SPACE]);
    end
end
