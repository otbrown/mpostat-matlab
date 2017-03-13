% EigenSolverTest.m
% Oliver Thomson Brown
% 2016-10-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true), matlab.unittest.fixtures.PathFixture('../../external/primme/Matlab')}) EigenSolverTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        effL = sparse([2, -1, -1, 1; -1i, 1i, -1, 1; 0, 1, 0, 2; 0, 1i, 0, 2i]);
    end

    methods (Test)
        function testThrowBadHermiticity(tc)
            tc.fatalAssertError(@()EigenSolver(tc.effL, true, 1E-9), ...
            'EigenSolver:badHermiticity');
        end

        function testHermitian(tc)
            HERMITIAN = true;
            H = ctranspose(tc.effL) * tc.effL;
            [vec, val] = EigenSolver(H, HERMITIAN, 1E-9);

            actualVal = 0;
            actualVec = [ 0.302322027807821 + 0.345148819050045i;
                          0.626057451236754 + 0.366562214671157i;
                         -0.334442121239489 + 0.140454316093354i;
                         -0.313028725618377 - 0.183281107335578i];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - vec), tc.absTol);
        end

        function testNonHermitian(tc)
            HERMITIAN = false;
            [vec, val] = EigenSolver(tc.effL, HERMITIAN, 1E-9);

            actualVal = 1.920620363439616 + 0.618180731783746i;
            actualVec = [0.909915536697958;
                         0.356041764971642;
                         0.150479197352233;
                         0.150479197352233];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - abs(vec)), tc.absTol);
        end
    end
end
