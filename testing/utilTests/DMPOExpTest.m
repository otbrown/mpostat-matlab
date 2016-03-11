% DMPOExpTest.m
% Oliver Thomson Brown
% 2016-03-11

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) DMPOExpTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        STATE;
        dmpo;
        prodDMPO;
        idOp;
        nOp;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3};
        testLENGTH = {7, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.STATE = tc.HILBY^tc.LENGTH - 1;
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.prodDMPO = ProdDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS, tc.STATE);
            tc.idOp = repmat(eye(tc.HILBY), [1, 1, tc.LENGTH]);
            n = zeros(tc.HILBY);
            n(tc.HILBY^2) = 1;
            tc.nOp = repmat(n, [1, 1, tc.LENGTH]);
        end
    end

    methods (Test)
        function testClass(tc)
            expect = DMPOExp(tc.dmpo, tc.idOp);
            tc.fatalAssertClass(expect, 'double');
        end

        function testArbTrace(tc)
            idExp = DMPOExp(tc.dmpo, tc.idOp);
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(idExp) - abs(tr));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testProdTrace(tc)
            idExp = DMPOExp(tc.prodDMPO, tc.idOp);
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(idExp) - abs(tr));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testNExpectation(tc)
            nExp = DMPOExp(tc.prodDMPO, tc.nOp);
            epsilon = abs(abs(nExp) - tc.LENGTH * (tc.HILBY - 1));
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end    
