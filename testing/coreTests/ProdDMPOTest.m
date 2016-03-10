% ProdDMPOTest.m
% Oliver Thomson Brown
% 2016-03-10

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) ProdDMPOTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        STATE;
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
            tc.STATE = randi([0, tc.HILBY^tc.LENGTH - 1]);
            tc.dmpo = ProdDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS, tc.STATE);
        end
    end

    methods (Test)
        % sufficient to test that trace = 1, and that the correct diagonal element
        % is 1 -- this means that all other combinations must be zero
        function testClass(tc)
            tc.fatalAssertClass(tc.dmpo, 'cell');
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
                colSz = size(tc.dmpo{site}, 2);
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

        function testStateElement(tc)
            local = FWBase(tc.STATE, tc.HILBY, tc.LENGTH) + 1;
            coefft = 1;
            for site = 1 : 1 : tc.LENGTH
                coefft = coefft * tc.dmpo{site}(:, :, local(site), local(site));
            end
            epsilon = abs(abs(coefft) - 1);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end
