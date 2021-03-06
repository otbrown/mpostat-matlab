% ZDMPOTest.m
% Oliver Thomson Brown
% 2016-03-10

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) ZDMPOTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        COMPRESS = 100;
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
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
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
            tc.fatalAssertError(@()ZDMPO(tc.HILBY, tc.LENGTH, BAD_COMPRESS), 'ZDMPO:BadCOMPRESS');
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
                colSz = size(tc.dmpo{site + 1}, 1);
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
            % since this normalisation happens after compression it *should* be
            % true, though other density matrix characteristics cannot be relied
            % upon
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(tr) - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end
