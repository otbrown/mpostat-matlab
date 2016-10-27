% DMPOTrace.m
% calculates the trace of a density matrix product operator
% Oliver Thomson Brown
% 2016-02-17
%
% [ trace ] = DMPOTrace(dmpo)
%
% RETURN
% trace:    (complex) double, the trace of the density matrix described
%           by dmpo -- ought to always be real-valued and one for
%           normalisation
%
% INPUT
% dmpo:     cell array, a density matrix product operator

function trace = DMPOTrace(dmpo)
    % pull constants
    LENGTH = size(dmpo, 1);
    HILBY = size(dmpo{1}, 3);

    % array allocate
    physCon = cell(LENGTH, 1);

    for site = 1 : 1 : LENGTH
        physCon{site} = zeros(size(dmpo{site}(:, :, 1, 1)));
        for state = 1 : 1 : HILBY
            physCon{site} = physCon{site} + dmpo{site}(:, :, state, state);
        end
    end

    trace = physCon{1};
    for site = 2 : 1 : LENGTH
        trace = trace * physCon{site};
    end
end
