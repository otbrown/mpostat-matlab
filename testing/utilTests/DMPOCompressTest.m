% DMPOCompressTest.m
% Oliver Thomson Brown
% 2016-02-26

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOCompressTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        mdmpo;
        pdmpo;
        cmdmpo;
        cpdmpo;
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
            % create density matrix mpo -- all at COMPRESS = 0
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.COMPRESS = testCOMPRESS;
            tc.mdmpo = MixDMPO(testHILBY, testLENGTH, 0);
            tc.pdmpo = ProdDMPO(testHILBY, testLENGTH, 0, 0);
            tc.COMPRESS = testCOMPRESS;
            tc.cmdmpo = DMPOCompress(tc.mdmpo, tc.COMPRESS, ...
                                        tc.HILBY, tc.LENGTH);
            tc.cmdmpo = TrNorm(tc.cmdmpo);
            tc.cpdmpo = DMPOCompress(tc.pdmpo, tc.HILBY^2, ...
                                        tc.HILBY, tc.LENGTH);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.cmdmpo, 'cell');
            tc.fatalAssertClass(tc.cpdmpo, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.cmdmpo, size(tc.mdmpo));
            tc.fatalAssertSize(tc.cpdmpo, size(tc.pdmpo));
        end

        function testTensorShape(tc)
            for site = 1 : 1 : tc.LENGTH
                mdim = min(size(tc.mdmpo{site}), tc.COMPRESS);
                pdim = min(size(tc.pdmpo{site}), tc.HILBY^2);

                tc.fatalAssertSize(tc.cmdmpo{site}, [mdim(1), mdim(2), ...
                    tc.HILBY, tc.HILBY]);
                tc.fatalAssertSize(tc.cpdmpo{site}, [pdim(1), pdim(2), ...
                    tc.HILBY, tc.HILBY]);
            end
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

                compCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    compCoefft = compCoefft ...
                                    * tc.cmdmpo{site}(:, :, bra, ket);
                end

                epsilon = abs((1/SPACE) - compCoefft);
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end

        function testProductTrace(tc)
            % the product state should be easily representable with
            % minimal dimensions so trace should be good in spite of
            % compression
            tr = DMPOTrace(tc.cpdmpo);
            epsilon = abs(tr - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testProductElement(tc)
            % again, 0s product state, so first element in DM should be 1
            coefft = 1;
            for site = 1 : 1 : tc.LENGTH
                coefft = coefft * tc.cpdmpo{site}(:, :, 1, 1);
            end
            epsilon = abs(coefft - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end
