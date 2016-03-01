% DMPOCompressTest.m
% Oliver Thomson Brown
% 2016-02-26

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) DMPOCompressTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        dmpo;
        cdmpo;
        cTr;
        prodState;
        systemSz;
        HILBY;
        COMPRESS;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 2, 3, 4};
        testSystemSz = {7, 4, 8, 5, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testSystemSz)
            % create density matrix mpo -- all at COMPRESS = 0
            tc.HILBY = testHILBY;
            tc.systemSz = testSystemSz;
            tc.prodState = randi([0, testHILBY^testSystemSz - 1]);
            tc.dmpo= ProdDMPO(testHILBY, testSystemSz, 0, tc.prodState);
            tc.COMPRESS = testHILBY^(testSystemSz - 2);
            tc.cdmpo = DMPOCompress(tc.dmpo, tc.COMPRESS);
            tc.cTr = DMPOTrace(tc.cdmpo);
        end
    end

    methods (Test)
        function testThrowBadCompress(tc)
            % have to anonymise DMPOCompress here to ensure that error is caught
            % by fatal assert
            BAD_COMPRESS = tc.HILBY^2 - 1;
            tc.fatalAssertError(@()DMPOCompress(tc.dmpo, BAD_COMPRESS), 'DMPOCompress:BadCOMPRESS');
        end

        function testClass(tc)
            tc.fatalAssertClass(tc.cdmpo, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.cdmpo, size(tc.dmpo));
        end

        function testTensorShape(tc)
            for site = 1 : 1 : tc.systemSz
                dims = min(size(tc.dmpo{site}), tc.COMPRESS);

                tc.fatalAssertSize(tc.cdmpo{site}, [dims(1), dims(2), tc.HILBY, tc.HILBY]);
            end
        end

        function testTraceUnity(tc)
            tc.assertLessThan(abs(abs(tc.cTr) - 1), tc.absTol);
        end

        function testTraceImag(tc)
            tc.assertLessThan(imag(tc.cTr), tc.absTol);
        end

        function testTraceDrift(tc)
            tr = DMPOTrace(tc.dmpo);
            drift = abs(abs(tr) - abs(tc.cTr));

            tc.assertLessThan(drift, tc.absTol);
        end
    end
end
