% ConvTestTest.m
% Oliver Thomson Brown
% 2016-09-08

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) ConvTestTest < matlab.unittest.TestCase

    properties
        THRESHOLD = 1E-7;
        LENGTH = 30;
        FALSE = false;
        TRUE = true;
    end

    methods (Test)
        function testTooFewUpdates(tc)
            eigTracker = NaN(tc.LENGTH, 1);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertTrue(isnan(convergence));
            tc.assertEqual(convFlag, tc.FALSE);

            eigTracker(1 : tc.LENGTH - 1) = 0;
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertTrue(isnan(convergence));
            tc.assertEqual(convFlag, tc.FALSE);
        end

        function testNotConverged(tc)
            eigTracker = rand(tc.LENGTH, 1);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertGreaterThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.FALSE);

            eigTracker = tc.LENGTH * tc.THRESHOLD  * eigTracker;
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertGreaterThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.FALSE);

            eigTracker(1 : 2 : tc.LENGTH) = -1 * eigTracker(1 : 2 : tc.LENGTH);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertGreaterThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.FALSE);
        end

        function testConverged(tc)
            eigTracker = ones(tc.LENGTH, 1);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);

            eigTracker = tc.THRESHOLD * rand(tc.LENGTH, 1);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);

            eigTracker(1 : 2 : tc.LENGTH) = -0.5 * eigTracker(1 : 2 : tc.LENGTH);
            [convFlag, convergence] = ConvTest(eigTracker, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);
        end
    end
end
