% GrowBlock.m
% interface function for GrowLeft and GrowRight
% Oliver Thomson Brown
% 2016-10-24
%
% updateBlock = GrowBlock(dmpo, mpo, left, right, site, direction)
%

function [updateBlock] = GrowBlock(dmpo, mpo, left, right, site, direction)
    % pull site tensors
    siteTensor = dmpo{site};
    mpoTensor = mpo{site};

    % pull dimensions
    [ROW_SIZE, COL_SIZE, HILBY, ~] = size(siteTensor);
    [~, ~, ~, ~, OP_ROW, OP_COL] = size(mpoTensor);

    % call correct Grow function
    if strcmpi(direction, 'L')
            lBlock = left{site};
            updateBlock = GrowLeft(siteTensor, mpoTensor, lBlock, ...
                                    ROW_SIZE, COL_SIZE, HILBY, OP_COL);
    elseif strcmpi(direction, 'R');
            rBlock = right{site};
            updateBlock = GrowRight(siteTensor, mpoTensor, rBlock, ...
                                    ROW_SIZE, COL_SIZE, HILBY, OP_ROW);
    else
        ME = MException('GrowBlock:BadDirection', 'GrowBlock accepts the case insensitive characters L or R to indicate direction. %s was entered.', direction);
        throw(ME);
    end
end
