% DMPOResizeTest.m
% Oliver Thomson Brown
% 2016-11-17

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOResizeTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        HILBY = 3;
        LENGTH = 6;
        OLD_COMPRESS = 80;
        ENLARGE;
        xdmpo;
        cdmpo;
        elDMPO;
        compDMPO;
    end

    methods (TestClassSetup)
        function ClassSetup(tc)
            tc.xdmpo = DDMPO(tc.HILBY, tc.LENGTH, 0);
            tc.cdmpo = DDMPO(tc.HILBY, tc.LENGTH, tc.OLD_COMPRESS);
            tc.ENLARGE = tc.OLD_COMPRESS + 20;
            tc.elDMPO = DMPOResize(tc.cdmpo, tc.ENLARGE);
            tc.compDMPO = DMPOResize(tc.xdmpo, tc.OLD_COMPRESS);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.elDMPO, 'cell');
            tc.fatalAssertClass(tc.compDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.elDMPO, size(tc.cdmpo));
            tc.fatalAssertSize(tc.compDMPO, size(tc.xdmpo));
        end

        function testThrowBadCompress(tc)
              % have to anonymise DMPOResize here to ensure that error is
              % caught by fatal assert
              COMPRESS = tc.HILBY^2 - 1;
              tc.fatalAssertError(@()DMPOResize(tc.xdmpo, ...
                  COMPRESS), 'DMPOResize:BadCOMPRESS');
        end

        function testTrace(tc)
            elTr = DMPOTrace(tc.elDMPO);
            compTr = DMPOTrace(tc.compDMPO);
            tc.assertLessThan(abs(elTr - 1), tc.absTol);
            tc.assertLessThan(abs(compTr - 1), tc.absTol);
        end

        function testSameCompress(tc)
            rsdmpo = DMPOResize(tc.cdmpo, tc.OLD_COMPRESS);
            tc.assertEqual(rsdmpo, tc.cdmpo);
        end

        function testEnlargeExact(tc)
            rsdmpo = DMPOResize(tc.xdmpo, tc.HILBY^tc.LENGTH);
            tc.assertEqual(rsdmpo, tc.xdmpo);
        end

        function testNotLarger(tc)
            for site = 1 : 1 : tc.LENGTH
                elDim = size(tc.elDMPO{site});
                compDim = size(tc.compDMPO{site});

                tc.assertLessThan(elDim, tc.ENLARGE+1);
                tc.assertLessThan(compDim, tc.OLD_COMPRESS+1)
            end
        end
    end
end
