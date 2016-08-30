% TrNormTest.m
% Oliver Thomson Brown
% 2016-03-09

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) TrNormTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        dmpo;
        normDMPO;
        tr;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5}
        testLENGTH = {7, 6, 5, 5}
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            % form a DMPO
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            % mess it up
            tc.dmpo = DMPOScalarDiv(tc.dmpo, 0.1);
            % trace norm it
            tc.normDMPO = TrNorm(tc.dmpo);
            % calculate new trace
            tc.tr = DMPOTrace(tc.normDMPO);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.normDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.normDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            for site = 1 : 1 : tc.LENGTH
                tc.fatalAssertSize(tc.normDMPO{site}, size(tc.dmpo{site}));
            end
        end

        function testTrace(tc)
            epsilon = abs(abs(tc.tr) - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end
