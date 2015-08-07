% SuperDMPSTest.m
% a testing class for SuperDMPS.m, which creates a density matrix MPS
% which represents a system in an even superposition state
% Oliver Thomson Brown
% 2015-08-07

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) SuperDMPSTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        superDMPS;
        systemSz;
        HILBY;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 2, 2, 3, 4, 5, 6};
        testSystemSz = {2, 3, 5, 7, 5, 4, 3, 3};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(testCase, testHILBY, testSystemSz)
            % create density matrix mps -- all at COMPRESS = 0
            testCase.superDMPS = SuperDMPS(testHILBY, testSystemSz, 0);
            testCase.HILBY = testHILBY;
            testCase.systemSz = testSystemSz;
        end
    end

    methods (Test)
        function testClass(testCase)
            testCase.assertClass(testCase.superDMPS, 'cell');
        end

        function testCellShape(testCase)
            testCase.assertSize(testCase.superDMPS, [testCase.systemSz, 1]);
        end

        function testMatrixShape(testCase)
            braKetSpace = (testCase.HILBY)^2;
            % initialise rowSz, colSz, dim3Sz
            rowSz = 1;
            colSz = braKetSpace;
            dim3Sz = braKetSpace;
            % sweep through each site in the MPS and assert the correct shape
            for site = 1 : 1 : testCase.systemSz
                testCase.assertSize(testCase.superDMPS{site}, [rowSz, colSz, dim3Sz]);
                % locate site in the system and determine correct size exponent
                lLen = site + 1;
                rLen = testCase.systemSz - site - 1;
                len = min(lLen, rLen);
                % set next site rowSz and colSz
                rowSz = colSz;
                colSz = testCase.HILBY^(2*len);
            end
        end

        function testCorrect(testCase)
            densityMatrix = DMRebuild(testCase.superDMPS);
            actualDM = zeros(testCase.HILBY^(2*testCase.systemSz), 1);
            actualDM(1 : (testCase.HILBY^testCase.systemSz + 1) : end) = 1 / (testCase.HILBY^testCase.systemSz);
            testCase.assertEqual(densityMatrix, actualDM, 'AbsTol', testCase.absTol);
        end
    end
end
