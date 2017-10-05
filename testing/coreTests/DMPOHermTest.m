% DMPOHermTest.m
% Oliver Thomson Brown
% 2016-03-09

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOHermTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        sample = 0.1;
        SAMPLE_MAX = 100;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        dmpo;
        hDMPO;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testLENGTH = {7, 6, 5, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            % build a DMPO
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            % (hopefully) make it Hermitian
            tc.hDMPO = DMPOHerm(tc.dmpo);
        end
    end

    methods (Test)
        % essentially the same set of tests here as DMPOConjTest, which returns
        % the Hermitian conjugate of a DMPO, only here we are checking the DMPO is
        % itself Hermitian
        function testClass(tc)
            tc.fatalAssertClass(tc.hDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.hDMPO, [tc.LENGTH, 1]);
        end

        function testTensorShape(tc)
            % since DMPOHerm involves adding the DMPO to itself, virtual
            % dimensions will double
            % first site
            dims = size(tc.dmpo{1});
            dims(2) = 2 * dims(2);
            tc.fatalAssertSize(tc.hDMPO{1}, dims);
            % between sites
            for site = 2 : 1 : tc.LENGTH - 1
                dims = size(tc.dmpo{site});
                dims(1 : 2) = 2 * dims(1 : 2);
                tc.fatalAssertSize(tc.hDMPO{site}, dims);
            end
            % last site
            dims = size(tc.dmpo{tc.LENGTH});
            dims(1) = 2 * dims(1);
            tc.fatalAssertSize(tc.hDMPO{tc.LENGTH}, dims);
        end

        function testTracePreservation(tc)
            tr1 = abs(DMPOTrace(tc.dmpo));
            tr2 = abs(DMPOTrace(tc.hDMPO));
            drift = abs(tr1 - tr2);
            tc.assertLessThan(drift, tc.absTol);
        end

        function testTraceReal(tc)
            trError = imag(DMPOTrace(tc.hDMPO));
            tc.assertLessThan(trError, tc.absTol);
        end

        function testHerm(tc)
            % samples the density matrix and checks each element's transpose is
            % it's conjugate
            SPACE = tc.HILBY^tc.LENGTH;
            sampleSz = min(floor(tc.sample * SPACE^2), tc.SAMPLE_MAX);
            for testNum = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                tCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = braBits(site);
                    coefft = coefft * tc.hDMPO{site}(:, :, bra, ket);
                    tCoefft = tCoefft * tc.hDMPO{site}(:, :, ket, bra);
                end

                epsilon = abs(coefft - conj(tCoefft));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
