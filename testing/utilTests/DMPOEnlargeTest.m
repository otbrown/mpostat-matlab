% DMPOEnlargeTest.m
% Oliver Thomson Brown
% 2016-11-16

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOEnlargeTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        mdmpo;
        pdmpo;
        bmdmpo;
        bpdmpo;
        LENGTH;
        HILBY;
        COMPRESS;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 2, 3};
        testLENGTH = {7, 4, 8, 5};
        testCOMPRESS = {48, 12, 128, 60};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH, testCOMPRESS)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.COMPRESS = testCOMPRESS;
            tc.mdmpo = MixDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.pdmpo = ProdDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.bmdmpo = DMPOEnlarge(tc.mdmpo, tc.COMPRESS + 10, ...
                                    tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.bpdmpo = DMPOEnlarge(tc.pdmpo, tc.COMPRESS + 10, ...
                                    tc.HILBY, tc.LENGTH, tc.COMPRESS);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.bmdmpo, 'cell');
            tc.fatalAssertClass(tc.bpdmpo, 'cell');
        end

        function testSystemSize(tc)
        end
    end
end
