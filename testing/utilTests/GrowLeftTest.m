% GrowLeftTest.m
% Oliver Thomson Brown
% 2016-03-22

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) GrowLeftTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        dmpo;
        leftBlock;
        impo;
        OP_ROW = 1;
        OP_COL = 1;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 3};
        testLENGTH = {3, 4, 3};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.dmpo = LCan(tc.dmpo, [1 : 1 : (tc.LENGTH - 1)]);
            tc.leftBlock = cell(tc.LENGTH, 1);
            tc.leftBlock{1} = 1;
            tc.impo = reshape(eye(tc.HILBY^2), [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            % populate leftBlock
            for site = 2 : 1 : tc.LENGTH
                [rowSz, colSz, ~, ~] = size(tc.dmpo{site - 1});
                tc.leftBlock{site} = GrowLeft(tc.dmpo{site - 1}, tc.impo, tc.leftBlock{site - 1}, rowSz, colSz, tc.HILBY, tc.OP_ROW, tc.OP_COL);
            end
        end
    end

    methods (Test)
        function vecNormTest(tc)
            % calculate vector norm by using LeftGrow past the end of the system
            [rowSz, colSz, ~, ~] = size(tc.dmpo{tc.LENGTH});
            testVecNorm = GrowLeft(tc.dmpo{tc.LENGTH}, tc.impo, tc.leftBlock{tc.LENGTH}, rowSz, colSz, tc.HILBY, tc.OP_ROW, tc.OP_COL);

            % calculate the vector norm the old fashioned way...
            % RETURN ALLOCATION
            SPACE = tc.HILBY^tc.LENGTH;
            rhoVec = zeros(SPACE^2, 1);

            % REBUILD VECTORISED DENSITY MATRIX
            for ketState = 0 : 1 : (SPACE - 1)
                for braState = 0 : 1 : (SPACE - 1)
                    stateDex = ketState * SPACE + braState + 1;
                    braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                    ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                    coefft = 1;
                    for site = 1 : 1 : tc.LENGTH
                        bra = braBits(site);
                        ket = ketBits(site);
                        coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    end

                    rhoVec(stateDex) = coefft;
                end
            end
            vecNorm = ctranspose(rhoVec) * rhoVec;

            epsilon = abs(testVecNorm - vecNorm);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end