% FWBaseTest.m
% Oliver Thomson Brown
% 2016-03-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) FWBaseTest < matlab.unittest.TestCase

    properties (TestParameter)
        n = {0, 1, 8, 17, 42, 100};
        BASE = {2, 2, 2, 3, 4, 10};
        WIDTH = {1, 1, 4, 4, 8, 3};
    end

    methods (Test)
        function testThrowNegN(tc)
            % have to anonymise function to ensure error is caught by test harness
            tc.fatalAssertError(@()FWBase(-10, 2, 4), 'FWBase:BadN');
        end

        function testThrowFloatN(tc)
            tc.fatalAssertError(@()FWBase(4.2, 2, 4), 'FWBase:BadN');
        end

        function testThrowBadWidth(tc)
            tc.fatalAssertError(@()FWBase(16, 2, 0), 'FWBase:BadWidth');
            tc.fatalAssertError(@()FWBase(16, 2, 2), 'FWBase:BadWidth');
        end
    end

    methods (Test, ParameterCombination='sequential')
        function testClass(tc, n, BASE, WIDTH)
            bits = FWBase(n, BASE, WIDTH);
            tc.fatalAssertClass(bits, 'double');
        end

        function testWidth(tc, n, BASE, WIDTH)
            bits = FWBase(n, BASE, WIDTH);
            tc.fatalAssertSize(bits, [WIDTH, 1]);
        end

        function testDecimal(tc, n, BASE, WIDTH)
            bits = FWBase(n, BASE, WIDTH);
            dec = 0;
            for bit = WIDTH : -1 : 1
                dec = dec + bits(bit) * BASE^(WIDTH - bit);
            end
            tc.assertEqual(dec, n);
        end
    end
end
