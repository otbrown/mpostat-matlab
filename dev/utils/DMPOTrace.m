% DMPOTrace.m
% calculates the trace of a supplied density matrix product operator
% Oliver Thomson Brown
% 2016-02-17
%
% trace = DMPOTrace(dmpo)
%
% RETURN
% trace     : complex double, the trace of the density matrix described by dmpo
%
% INPUTS
% dmpo      : cell array, an arbitrary density matrix product operator

function trace = DMPOTrace(dmpo)
    % we just want the coefficients that correspond to both input and output
    % states of the density matrix being the same (the diagonal elements)

    % gather data
    LENGTH = size(dmpo, 1);
    HILBY = size(dmpo{1}, 3);
    SPACE = HILBY^LENGTH - 1;

    % initialise return
    trace = 0;

    % loop through states in the system, find their coefft and add it in
    for state = 0 : 1 : SPACE
        bits = FWBase(state, HILBY, LENGTH) + 1;
        coefft = 1;
        for site = 1 : 1 : LENGTH
            coefft = coefft * dmpo{site}(:, :, bits(site), bits(site));
        end
        trace = trace + coefft;
    end
end
