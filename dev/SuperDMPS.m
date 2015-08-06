% SuperDMPS.m
% function which forms a density matrix MPO for a system in an even super-position of it's principle states
% Oliver Thomson Brown
% 2015-08-06
%
% Return:
% densityMPS    : cell array, L * 1, contains the matrix product state, physical indices in the third dimension 
%               -- each site has a ket and a bra, ket changes first
%
% Inputs:
% HILBY     : int, dimension of the local state-space 
% LENGTH    : int, number of the sites in the chaini
% COMPRESS  : int, maximum size of on-site matrices, supplying 0 returns an exact MPO

function [densityMPS] = SuperDMPS(HILBY, LENGTH, COMPRESS)
    if COMPRESS == 0
        COMPRESS = Inf;
    end

    % RETURN ALLOCATION
    densityMPS = cell(LENGTH, 1);

    % FIRST AND LAST SITES
    densityMPS{1} = zeros(1, HILBY^2, HILBY^2);
    densityMPS{LENGTH} = zeros(HILBY^2, 1, HILBY^2);

    trNorm = sqrt(HILBY^LENGTH);    
    trMap = permute(kron(ones(HILBY,1), eye(1, HILBY^2)), [3,2,1]);

    densityMPS{1}(:, :, 1 : (HILBY+1) : end) = trMap / trNorm;
    densityMPS{LENGTH}(:, :, 1 : (HILBY+1) : end) = trMap / trNorm;

    rowSize = HILBY^2;
    if LENGTH == 3
        colSize = HILBY^2;
    else 
        colSize = min(HILBY^4, COMPRESS);
    end

    for site = 2 : 1 : LENGTH - 1
        lLen = site + 1;
        rLen = LENGTH - site - 1;
        len = min(lLen, rLen);

        densityMPS{site} = ones(rowSize, colSize, HILBY^2);

        rowSize = colSize;
        colSize = min(COMPRESS, HILBY^(2*len));
    end
end
