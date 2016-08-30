% DMPOConjTest.m
% Oliver Thomson Brown
% 2016-03-01

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOConjTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        sample = 0.1;
        SAMPLE_MAX = 100;
        dmpo;
        conjDMPO;
        HILBY;
        systemSz;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testSystemSz = {7, 5, 4, 4};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testSystemSz)
            tc.HILBY = testHILBY;
            tc.systemSz = testSystemSz;
            tc.dmpo = DMPO(testHILBY, testSystemSz, 0);
            tc.conjDMPO = DMPOConj(tc.dmpo);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.conjDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.conjDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            for site = 1 : 1 : tc.systemSz
                tc.fatalAssertSize(tc.conjDMPO{site}, size(tc.dmpo{site}));
            end
        end

        function testTracePreservation(tc)
            tr = DMPOTrace(tc.dmpo);
            conjTr = DMPOTrace(tc.conjDMPO);
            drift = abs(tr - conjTr);
            tc.assertLessThan(drift, tc.absTol);
        end

        function testConj(tc)
            % samples the density matrix and checks it is conjugated
            SPACE = tc.HILBY^tc.systemSz;
            sampleSz = min(floor(tc.sample * (SPACE^2)), tc.SAMPLE_MAX);
            for testNum = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.systemSz) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.systemSz) + 1;

                coefft = 1;
                conjCoefft = 1;
                for site = 1 : 1 : tc.systemSz
                    bra = braBits(site);
                    ket = braBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    conjCoefft = conjCoefft * tc.conjDMPO{site}(:, :, ket, bra);
                end

                epsilon = abs(coefft - conj(conjCoefft));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
