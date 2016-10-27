% LCanTest.m
% Oliver Thomson Brown
% 2016-03-15

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) LCanTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        HILBY;
        LENGTH;
        COMPRESS = 50;
        SAMPLE_MAX = 100;
        dmpo;
        canDMPO;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testLENGTH = {7, 5, 4, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.canDMPO = tc.dmpo;
            for site = 1 : 1 : (tc.LENGTH - 1)
                siteTensor = tc.canDMPO{site};
                nextSiteTensor = tc.canDMPO{site+1};
                [ROW_SIZE, COL_SIZE, ~, ~] = size(siteTensor);
                NEXT_COL = size(nextSiteTensor, 2);
                [tc.canDMPO{site}, tc.canDMPO{site+1}] = LCan(siteTensor, ...
                    nextSiteTensor, tc.HILBY, ROW_SIZE, COL_SIZE, NEXT_COL);
            end
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.canDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.canDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            % LCan may reduce virtual dimensions, but should still be consistent
            rowSz = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                colSz = size(tc.canDMPO{site + 1}, 1);
                tc.fatalAssertSize(tc.canDMPO{site}, [rowSz, colSz, tc.HILBY, tc.HILBY]);
                rowSz = colSz;
            end
        end

        function testTracePreservation(tc)
            tr1 = DMPOTrace(tc.dmpo);
            tr2 = DMPOTrace(tc.canDMPO);
            epsilon = abs(abs(tr1) - abs(tr2));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testCan(tc)
            % should be able to re-arrange each tensor back to the unitary from
            % which it was formed -- A' * A = I
            for site = 1 : 1 : tc.LENGTH - 1
                [rowSz, colSz, ~, ~] = size(tc.canDMPO{site});
                M = reshape(tc.canDMPO{site}, [rowSz, colSz, tc.HILBY^2]);
                M = permute(M, [1, 3, 2]);
                M = reshape(M, [rowSz * tc.HILBY^2, colSz]);

                I = eye(colSz, colSz);
                epsilon = (ctranspose(M) * M) - I;
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end

        function testElements(tc)
            % LCan should not alter the density matrix being represented
            SPACE = tc.HILBY^tc.LENGTH;
            sampleSz = min(floor(0.1 * SPACE^2), tc.SAMPLE_MAX);
            for testNum = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                canCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    canCoefft = canCoefft * tc.canDMPO{site}(:, :, bra, ket);
                end

                epsilon = abs(abs(coefft) - abs(canCoefft));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
