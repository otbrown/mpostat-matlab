% DMPOScalarDivTest.m
% Oliver Thomson Brown
% 2016-03-02

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) DMPOScalarDivTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        sampleSz = 200;
        HILBY = 3;
        LENGTH = 7;
        COMPRESS = 100;
        scalar = 10;
        dmpo;
        divDMPO;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.divDMPO = DMPOScalarDiv(tc.dmpo, tc.scalar);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.divDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.divDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            for site = 1 : 1 : tc.LENGTH
                tc.fatalAssertSize(tc.divDMPO{site}, size(tc.dmpo{site}));
            end
        end

        function testTraceDiv(tc)
            tr = DMPOTrace(tc.dmpo);
            divTr = DMPOTrace(tc.divDMPO);
            epsilon = abs(abs(divTr) - abs(tr / tc.scalar));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testDiv(tc)
            SPACE = tc.HILBY^tc.LENGTH;
            for testNum = 1 : 1 : tc.sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                divCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    divCoefft = divCoefft * tc.divDMPO{site}(:, :, bra, ket);
                end

                epsilon = abs(abs(divCoefft) - abs(coefft / tc.scalar));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
