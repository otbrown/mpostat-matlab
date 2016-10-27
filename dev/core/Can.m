% Can.m
% interface function for left and right canonisation functions
% Oliver Thomson Brown
% 2016-10-24
%
% [ cdmpo ] = Can(dmpo, route, direction)
%
% RETURN
% cdmpo:    cell, a density matrix product operator with specified
%           tensor(s) brought into left- or right-canonical form
%
% INPUT
% dmpo:         cell, a density matrix product operator
% route:        double, the site or set of sites which should be brought
%               into canonical form
% direction:    character, 'L' or 'R', the direction of the canonisation
%               procedure

function [ cdmpo ] = Can(dmpo, route, direction)
    LENGTH = length(dmpo);
    HILBY = size(dmpo{1}, 3);

    cdmpo = dmpo;

    if strcmpi(direction, 'L')
        % check route is not invalid in a way Matlab won't notice
        if route(end) >= LENGTH
            msgID = 'Can:LCan:BadRoute';
            msg = sprintf(['Route cannnot extend to (or exceed) the ', ...
                            'last site in the system. System has %d ', ...
                            'sites, your route ended at %d.'], ...
                            LENGTH, route(end));
            ME = MException(msgID, msg);
            throw(ME);
        elseif route(1) > route(end)
            msgID = 'Can:LCan:BadRoute';
            msg = sprintf(['This route appears to be for the wrong ', ...
                            'direction. Left canonisation must be ', ...
                            'performed moving from the first site to ', ...
                            'the last.']);
            ME = MException(msgID, msg);
            throw(ME);
        end

        % perform left canonisation
        for site = route
            siteTensor = cdmpo{site};
            nextSiteTensor = cdmpo{site + 1};
            [ROW_SIZE, COL_SIZE, ~, ~] = size(siteTensor);
            NEXT_COL = size(nextSiteTensor, 2);
            [cdmpo{site}, cdmpo{site+1}] =  ...
                                    LCan(siteTensor, nextSiteTensor, ...
                                    HILBY, ROW_SIZE, COL_SIZE, NEXT_COL);
        end
    elseif strcmpi(direction, 'R')
        % check route is not invalid in a way Matlab won't notice
        if route(end) <= 1
            msgID = 'Can:RCan:BadRoute';
            msg = sprintf(['Route cannnot extend to (or exceed) the ', ...
                            'first site in the system. Your route ', ...
                            'ended at %d.'], route(end));
            ME = MException(msgID, msg);
            throw(ME);
        elseif route(1) < route(end)
            msgID = 'Can:RCan:BadRoute';
            msg = sprintf(['This route appears to be for the wrong ', ...
                            'direction. Right canonisation must be ', ...
                            'performed moving from the last site to ', ...
                            'the first.']);
            ME = MException(msgID, msg);
            throw(ME);
        end

        % perform right canonisation
        for site = route
            siteTensor = cdmpo{site};
            nextSiteTensor = cdmpo{site - 1};
            [ROW_SIZE, COL_SIZE, ~, ~] = size(siteTensor);
            NEXT_ROW = size(nextSiteTensor, 1);
            [cdmpo{site}, cdmpo{site-1}] = ...
                                    RCan(siteTensor, nextSiteTensor, ...
                                    HILBY, ROW_SIZE, COL_SIZE, NEXT_ROW);
        end
    else
        msgID = 'Can:BadDirection';
        msg = sprintf(['Direction should be case insensitive ', ...
                        'character L or R. Supplied direction was %s'], ...
                        direction);
        ME = MException(msgID, msg);
        throw(ME);
    end
end
