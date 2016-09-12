% EffOpTensorTest.m
% Oliver Thomson Brown
% 2016-05-04

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) EffOpTensorTest < matlab.unittest.TestCase

    properties
        HILBY;
        LENGTH;
        COMPRESS = 0;
        dmpo;
        left;
        impo;
        right;
        TEST_SITE;
        effTen;
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
            tc.impo = reshape(eye(tc.HILBY^2), [tc.HILBY, tc.HILBY, ...
                        tc.HILBY, tc.HILBY]);
            tc.left = cell(tc.LENGTH, 1);
            tc.right = cell(tc.LENGTH, 1);
            tc.left{1} = 1;
            tc.right{tc.LENGTH} = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{site});
                tc.left{site + 1} = GrowLeft(tc.dmpo{site}, tc.impo, ...
                tc.left{site}, ROW_SIZE, COL_SIZE, tc.HILBY, 1, 1);
            end
            for site = tc.LENGTH : -1 : 2
                [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{site});
                tc.right{site - 1} = GrowRight(tc.dmpo{site}, tc.impo, ...
                tc.right{site}, ROW_SIZE, COL_SIZE, tc.HILBY, 1, 1);
            end

            [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{tc.TEST_SITE});
            tc.effTen = EffOpTensor(tc.left{tc.TEST_SITE}, tc.impo, ...
            tc.right{tc.TEST_SITE}, ROW_SIZE, COL_SIZE, tc.HILBY);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.effTen, 'double');
        end

        function testTensorShape(tc)
            [ROW_SIZE, COL_SIZE, ~, ~] = size(tc.dmpo{tc.TEST_SITE});
            tc.fatalAssertSize(tc.effTen, [ROW_SIZE, COL_SIZE, COL_SIZE, ...
            ROW_SIZE, tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
        end

        function testNumNonZeroes(tc)
            numNonZeroes = nnz(tc.effTen);
            tc.assertEqual(numNonZeroes, tc.HILBY^2);
        end

        function testOnes(tc)
            for bra = 1 : 1 : tc.HILBY
                for ket = 1 : 1 : tc.HILBY
                    tc.assertEqual(tc.effTen(1,1,1,1,bra,ket,bra,ket), 1);
                end
            end
        end
    end
end
