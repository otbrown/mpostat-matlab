% GrowRight.m
% Oliver Thomson Brown
% 2016-03-22

function [updateBlock] = GrowRight(siteTensor, mpo, rightBlock, ROW_SIZE, COL_SIZE, HILBY, OP_ROW, OP_COL)
    updateBlock = zeros(ROW_SIZE, OP_ROW, ROW_SIZE);

    conjTensor = zeros(COL_SIZE, ROW_SIZE, HILBY, HILBY);
    for bra = 1 : 1 : HILBY
        for ket = 1 : 1 : HILBY
            conjTensor(:, :, bra, ket) = ctranspose(siteTensor(:, :, bra, ket));
        end
    end

    for conjCol = 1 : 1 : ROW_SIZE
        for opRow = 1 : 1 : OP_ROW
            for row = 1 : 1 : ROW_SIZE
                BWFB = 0;
                for conjBra = 1 : 1 : HILBY
                    for conjKet = 1 : 1 : HILBY
                        for conjRow = 1 : 1 : COL_SIZE
                            WFB = 0;
                            for bra = 1 : 1 : HILBY
                                for ket = 1 : 1 : HILBY
                                    for opCol = 1 : 1 : OP_COL
                                        FB = 0;
                                        for col = 1 : 1 : COL_SIZE
                                            FB = FB + rightBlock(conjRow, opCol, col) * siteTensor(row, col, bra, ket);
                                        end
                                        WFB = WFB + mpo(bra, ket, conjBra, conjKet, opRow, opCol) * FB;
                                    end
                                end
                            end
                            BWFB = BWFB + conjTensor(conjRow, conjCol, conjBra, conjKet) * WFB;
                        end
                    end
                end
                updateBlock(conjCol, opRow, row) = BWFB;
            end
        end
    end 
end
