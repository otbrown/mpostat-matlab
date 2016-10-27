
% EffLTest.m
% Oliver Thomson Brown
% 2016-09-29

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) EffLTest < matlab.unittest.TestCase

    properties
        HILBY = 2;
        LENGTH = 5;
        COMPRESS = 0;
        LDIM;
        effL;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            TARGET = ceil(tc.LENGTH / 2);
            dmpo = MixDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            impo = cell(tc.LENGTH, 1);
            for site = 1 : 1 : tc.LENGTH
                impo{site} = reshape(eye(tc.HILBY^2), ...
                                    [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            end
            left = cell(tc.LENGTH, 1);
            left{1} = 1;
            right = cell(tc.LENGTH, 1);
            right{tc.LENGTH} = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpo{site});
                left{site + 1} = GrowLeft(dmpo{site}, ...
                impo{site}, left{site}, ROW_SIZE, COL_SIZE, ...
                tc.HILBY, 1);
            end
            for site = tc.LENGTH : -1 : 2
                [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpo{site});
                right{site - 1} = GrowRight(dmpo{site}, ...
                impo{site}, right{site}, ROW_SIZE, COL_SIZE, ...
                tc.HILBY, 1);
            end

            [ROW_SIZE, COL_SIZE, ~, ~] = size(dmpo{TARGET});
            tc.LDIM = ROW_SIZE * COL_SIZE * tc.HILBY^2;

            tc.effL = EffL(TARGET, dmpo, impo, left, right);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.effL, 'double');
        end

        function testSparseReturn(tc)
            tc.assertTrue(issparse(tc.effL));
        end

        function testSize(tc)
            tc.assertSize(tc.effL, [tc.LDIM, tc.LDIM]);
        end
    end
end
