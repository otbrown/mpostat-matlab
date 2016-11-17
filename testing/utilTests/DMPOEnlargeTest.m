% DMPOEnlargeTest.m
% Oliver Thomson Brown
% 2016-11-16

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOEnlargeTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        dmpo;
        bdmpo;
        LENGTH;
        HILBY;
        COMPRESS;
        NEW_COMPRESS;
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
            tc.NEW_COMPRESS = testCOMPRESS + 10;
            tc.dmpo = MixDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.bdmpo = DMPOEnlarge(tc.dmpo, tc.NEW_COMPRESS, ...
                                    tc.HILBY, tc.LENGTH);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.bdmpo, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.bdmpo, size(tc.dmpo));
        end

        function testTensorShape(tc)
            rowSz = 1;
            for site = 1 : 1 : tc.LENGTH
                if site < ceil(tc.LENGTH/2)
                    len = site;
                else
                    len = tc.LENGTH - site;
                end
                colSz = min(tc.HILBY^(2*len), tc.NEW_COMPRESS);

                tc.fatalAssertSize(tc.bdmpo{site}, [rowSz, colSz, ...
                                    tc.HILBY, tc.HILBY]);
                rowSz = colSz;
            end
        end

        function testTrace(tc)
            tr = DMPOTrace(tc.bdmpo);
            epsilon = abs(tr - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testMixedElements(tc)
            % trace normalisation means every element should be 1/SPACE
            % -- the density matrix is just proportional to a ones matrix
            SPACE = tc.HILBY^tc.LENGTH;
            sampleSz = min(floor(0.1 * SPACE^2), 100);
            for test = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.bdmpo{site}(:, :, bra, ket);
                end

                epsilon = abs((1/SPACE) - coefft);
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
