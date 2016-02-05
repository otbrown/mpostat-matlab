% SuperDMPS.m
% function which returns a DMPS which represents
% an even superposition of states -- not practically
% useful, but helpful for checking chosen structure
% of the DMPS is sound, and that in principle one can
% recover the original density matrix
% Oliver Thomson Brown
% 2016-02-01
%
% densityMPO = SuperDMPO(HILBY, LENGTH, COMPRESS)
%
% [RETURN]
% densityMPO:   LENGTH x 1 cell array, format is densityMPO{site}(row, column, braState, ketState)
%				as in < braState | rho | ketState>
%
% [INPUTS]
% HILBY:        integer, the dimension of the local state space (i.e. for qubits, HILBY = 2)
% LENGTH:       integer, the size of the 1-D system
% COMPRESS:     integer, the maximum dimension of the matrices, enter 0 for an uncompressed MPS

function [densityMPO] = SuperDMPO(HILBY, LENGTH, COMPRESS)

    % COMPRESS == 0 means no compression
    if COMPRESS == 0
        COMPRESS = Inf;
    elseif COMPRESS < HILBY^2
        fprintf(1, 'Maximum matrix size cannot be less than the size of the 2-site density matrix. Please set COMPRESS >= %d.\n', HILBY^2);
        densityMPO = 0;
        return;
    end 

    % return allocation
    densityMPO = cell(LENGTH, 1);

    % first and last site
    norm = 1 / sqrt(HILBY^LENGTH);
    
    densityMPO{1} = zeros(1, HILBY^2, HILBY, HILBY);
    densityMPO{LENGTH} = zeros(HILBY^2, 1, HILBY, HILBY);

    densityMPO{1}(:,1,:,:) = norm;
    densityMPO{LENGTH}(:,1,:,:) = norm;

    colSz = HILBY^2;

    for site  = 2 : 1 : LENGTH - 1
        if site < ceil(LENGTH / 2)
            len = site;
        else
            len = LENGTH - site;
        end

        rowSz = colSz;
        colSz = min(HILBY^(2*len), COMPRESS);

        densityMPO{site} = zeros(rowSz, colSz, HILBY, HILBY);
        for i = 1 : 1 : min(rowSz, colSz)
            densityMPO{site}(i, i, :, :) = 1;
        end

    end
end
