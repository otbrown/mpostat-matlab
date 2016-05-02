% SuperDMPOTest.m
% Oliver Thomson Brown
% 2016-03-14

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) SuperDMPOTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        COMPRESS = 100;
        SAMPLE_MAX = 100;
        HILBY;
        LENGTH;
        dmpo;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testLENGTH = {7, 5, 4, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.dmpo = SuperDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.dmpo, 'cell');
        end

        function testThrowBadCompress(tc)
            % have to anonymise function to ensure that error is caught
            % by fatal assert
            BAD_COMPRESS = tc.HILBY^2 - 1;
            tc.fatalAssertError(@()MDMPO(tc.HILBY, tc.LENGTH, BAD_COMPRESS), 'MDMPO:BadCOMPRESS');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.dmpo, [tc.LENGTH, 1]);
        end

        function testTensorShape(tc)
            % since this function forms the DMPO virtual dimensions are tested for
            % validity, rather than 'correctness' -- that is that each tensor can
            % be multiplied into the next
            rowSz = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                colSz = size(tc.dmpo{site+1}, 1);
                tc.fatalAssertSize(tc.dmpo{site}, [rowSz, colSz, tc.HILBY, tc.HILBY]);
                % update rowSz for next site
                rowSz = colSz;
            end
            tc.fatalAssertSize(tc.dmpo{tc.LENGTH}, [rowSz, 1, tc.HILBY, tc.HILBY]);
        end

        function testCompression(tc)
            for site = 1 : 1 : tc.LENGTH
                tc.assertLessThan(size(tc.dmpo{site}), tc.COMPRESS + 1);
            end
        end

        function testTraceUnity(tc)
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(tr) - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testElements(tc)
            % similar to the method used to test DMPOHerm, and DMPOConj, the
            % density matrix is sampled, and it is checked that the elements equal
            % 1 / (HILBY^LENGTH)
            SPACE = tc.HILBY^tc.LENGTH;
            sampleSz = min(floor(0.1 * SPACE^2), tc.SAMPLE_MAX);
            for testNum = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                end

                epsilon = abs(abs(coefft) - (1 / SPACE));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
