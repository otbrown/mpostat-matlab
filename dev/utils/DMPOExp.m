% DMPOExp.m
% function which calculates an expectation value from an operator and a
% density matrix product operator <O> = Tr{O*rho}
% Oliver Thomson Brown
% 2016-03-11
%
% [ expect ] = DMPOExp(dmpo, op)
%
% RETURN
% expect:   (complex) double, the expectation value Tr{O*rho}
%
% INPUT
% dmpo:     cell array, a density matrix product operator
% op:       double array, op(ketState, braState, site);

function [expect] = DMPOExp(dmpo, op)
    % pull constants
    LENGTH = size(dmpo, 1);
    HILBY = size(dmpo{1}, 3);

    % array allocate
    physCon = cell(LENGTH, 1);
    for site = 1 : 1 : LENGTH
        physCon{site} = zeros(size(dmpo{site}(:, :, 1, 1)));
        for bra = 1 : 1 : HILBY
            for ket = 1 : 1 : HILBY
                physCon{site} = physCon{site} + op(ket, bra, site) * dmpo{site}(:, :, bra, ket);
            end
        end
    end

    expect = physCon{1};
    for site = 2 : 1 : LENGTH
        expect = expect * physCon{site};
    end
end
