% DMPOSumTest.m
% Oliver Thomson Brown
% 2016-03-02

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOSumTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        sampleSz = 100;
        HILBY = 3;
        LENGTH = 8;
        COMPRESS = 100;
        dmpo;
        sumDMPO;
        badHilbyDMPO;
        badLengthDMPO;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.sumDMPO = DMPOSum(tc.dmpo, tc.dmpo);
            tc.badHilbyDMPO = DMPO(tc.HILBY - 1, tc.LENGTH, tc.COMPRESS);
            tc.badLengthDMPO = DMPO(tc.HILBY, tc.LENGTH - 1, tc.COMPRESS);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.sumDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.sumDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            % virtual dimensions should have doubled under addition
            for site = 2 : 1 : tc.LENGTH - 1
                dims = size(tc.dmpo{site});
                dims(1:2) = 2 * dims(1:2);
                tc.fatalAssertSize(tc.sumDMPO{site}, dims);
            end

            % as ever first and last site are special snowflakes
            % first site
            dims = size(tc.dmpo{1});
            dims(2) = 2 * dims(2);
            tc.fatalAssertSize(tc.sumDMPO{1}, dims);
            % last site
            dims = size(tc.dmpo{tc.LENGTH});
            dims(1) = 2 * dims(1);
            tc.fatalAssertSize(tc.sumDMPO{tc.LENGTH}, dims);
        end

        function testThrowBadHilby(tc)
            % have to anonymise function to ensure error is caught by test harness
            tc.fatalAssertError(@()DMPOSum(tc.dmpo, tc.badHilbyDMPO), 'DMPOSum:BadArgs');
        end

        function testThrowBadLength(tc)
            % have to anonymise function to ensure error is caught by test harness
            tc.fatalAssertError(@()DMPOSum(tc.dmpo, tc.badLengthDMPO), 'DMPOSum:BadArgs');
        end

        function testTraceSum(tc)
            tr = DMPOTrace(tc.dmpo);
            sumTr = DMPOTrace(tc.sumDMPO);
            epsilon = abs(abs(sumTr) - 2 * abs(tr));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testSum(tc)
            SPACE = tc.HILBY^tc.LENGTH;
            for testNum = 1 : 1 : tc.sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                sumCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    sumCoefft = sumCoefft * tc.sumDMPO{site}(:, :, bra, ket);
                end

                epsilon = abs(abs(sumCoefft) - 2 * abs(coefft));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
