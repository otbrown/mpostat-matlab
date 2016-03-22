% GrowLeft.m
% Oliver Thomson Brown
% 2016-03-22

function [updateBlock] = GrowLeft(siteTensor, mpo, leftBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)

updateBlock = zeros(COL_SIZE, OP_COL, COL_SIZE);
conjTensor = zeros(COL_SIZE, ROW_SIZE, HILBY, HILBY);
for bra = 1 : 1 : HILBY
    for ket = 1 : 1 : HILBY
        conjTensor(:, :, bra, ket) = ctranspose(siteTensor(:, :, bra, ket));
    end
end

for conjRow = 1 : 1 : COL_SIZE
    for opCol = 1 : 1 : OP_COL
        for col = 1 : 1 : COL_SIZE
            AWFA = 0;
            for conjBra = 1 : 1 : HILBY
                for conjKet = 1 : 1 : HILBY
                    for conjCol = 1 : 1 : ROW_SIZE
                        WFA = 0;
                        for bra = 1 : 1 : HILBY
                            for ket = 1 : 1 : HILBY
                                for opRow = 1 : 1 : OP_ROW
                                    FA = 0;
                                    for row = 1 : 1 : ROW_SIZE
                                        FA = FA + leftBlock(conjCol, opRow, row) * siteTensor(row, col, bra, ket);
                                    end
                                    WFA = WFA + mpo(bra, ket, conjBra, conjKet, opRow, opCol) * FA;
                                end
                            end
                        end
                        AWFA = AWFA + conjTensor(conjRow, conjCol, conjBra, conjKet) * WFA;
                    end
                end
            end
            updateBlock(conjRow, opCol, col) = AWFA;
        end
    end
end

end
