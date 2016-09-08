% ConvTestTest.m
% Oliver Thomson Brown
% 2016-09-08

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) ConvTestTest < matlab.unittest.TestCase

    properties
        THRESHOLD = 1E-7;
        eigTracker = NaN(100, 1);
        sampleSize = 25;
        FALSE = 0;
        TRUE = 1;
    end

    methods (Test)
        function testTooFewUpdates(tc)
            [convFlag, convergence] = ConvTest(tc.eigTracker, tc.sampleSize - 1, tc.sampleSize, tc.THRESHOLD);
            tc.assertTrue(isnan(convergence));
            tc.assertEqual(convFlag, tc.FALSE);

            [convFlag, convergence] = ConvTest(tc.eigTracker, tc.sampleSize, tc.sampleSize, tc.THRESHOLD);
            tc.assertTrue(isnan(convergence));
            tc.assertEqual(convFlag, tc.FALSE);
        end

        function testNotConverged(tc)
            eigTracker = tc.eigTracker;
            updCount = 2 * tc.sampleSize;
            eigTracker(1 : updCount) = rand(updCount, 1);
            [convFlag, convergence] = ConvTest(eigTracker, updCount, tc.sampleSize, tc.THRESHOLD);
            tc.assertGreaterThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.FALSE);

            eigTracker(1 : updCount) = ones(updCount, 1);
            eigTracker(updCount) = 1 + tc.THRESHOLD;
            [convFlag, convergence] = ConvTest(eigTracker, updCount, tc.sampleSize, tc.THRESHOLD);
            tc.assertGreaterThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.FALSE);
        end

        function testConverged(tc)
            eigTracker = tc.eigTracker;
            updCount = 2 * tc.sampleSize;
            eigTracker(1 : updCount) = ones(updCount, 1);
            [convFlag, convergence] = ConvTest(eigTracker, updCount, tc.sampleSize, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);

            eigTracker(1 : updCount) = ones(updCount, 1) + tc.THRESHOLD * rand(updCount, 1);
            [convFlag, convergence] = ConvTest(eigTracker, updCount, tc.sampleSize, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);

            eigTracker(1 : updCount) = zeros(updCount, 1) + (tc.THRESHOLD/10) * rand(updCount, 1);
            eigTracker(1 : 2 : updCount) = -1 * eigTracker(1 : 2 : updCount);
            [convFlag, convergence] = ConvTest(eigTracker, updCount, tc.sampleSize, tc.THRESHOLD);
            tc.assertLessThan(convergence, tc.THRESHOLD);
            tc.assertEqual(convFlag, tc.TRUE);
        end
    end
end
