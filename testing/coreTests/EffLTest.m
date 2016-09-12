% EffLTest.m
% Oliver Thomson Brown
% 2016-05-05

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) EffLTest < matlab.unittest.TestCase

    properties
        HILBY;
        LENGTH;
        COMPRESS = 0;
        dmpo;
        left;
        impo;
        right;
        TEST_SITE;
        effL;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3}
        testLENGTH = {3, 4}
    end

    methods (TestMethodSetup)
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.TEST_SITE = ceil(tc.LENGTH / 2);
            tc.dmpo = ProdDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS, 0);
            tc.impo = cell(tc.LENGTH, 1);
            for site = 1 : 1 : tc.LENGTH
                tc.impo{site} = reshape(eye(tc.HILBY^2), ...
                                [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            end
            tc.left = cell(tc.LENGTH, 1);
            tc.left{1} = 1;
            tc.right = cell(tc.LENGTH, 1);
            tc.right{tc.LENGTH} = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{site});
                tc.left{site + 1} = GrowLeft(tc.dmpo{site}, ...
                tc.impo{tc.TEST_SITE}, tc.left{site}, ROW_SIZE, COL_SIZE, ...
                tc.HILBY, 1, 1);
            end
            for site = tc.LENGTH : -1 : 2
                [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{site});
                tc.right{site - 1} = GrowRight(tc.dmpo{site}, ...
                tc.impo{tc.TEST_SITE}, tc.right{site}, ROW_SIZE, COL_SIZE, ...
                tc.HILBY, 1);
            end

            [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{tc.TEST_SITE});
            tc.effL = EffL(tc.TEST_SITE, tc.dmpo, tc.impo, tc.left, tc.right);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.effL, 'double');
        end

        function testShape(tc)
            [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{tc.TEST_SITE});
            tc.fatalAssertSize(tc.effL, [ROW_SIZE*COL_SIZE*tc.HILBY*tc.HILBY,...
                                    ROW_SIZE*COL_SIZE*tc.HILBY*tc.HILBY]);
        end

        function testNumNonZeroes(tc)
            numNonZeroes = nnz(tc.effL);
            tc.assertEqual(numNonZeroes, tc.HILBY^2);
        end

        function testOnes(tc)
            [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{tc.TEST_SITE});
            for bra = 0 : 1 : (tc.HILBY - 1)
                for ket = 0 : 1 : (tc.HILBY - 1)
                    joindex = ket * tc.HILBY * ROW_SIZE * COL_SIZE ...
                            + bra * ROW_SIZE * COL_SIZE + 1;
                    tc.assertEqual(tc.effL(joindex, joindex), 1);
                end
            end
        end
    end
end
