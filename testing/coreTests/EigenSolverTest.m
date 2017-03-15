% EigenSolverTest.m
% Oliver Thomson Brown
% 2016-10-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true), matlab.unittest.fixtures.PathFixture('../../external/primme/Matlab')}) EigenSolverTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        effL = sparse([2, -1, -1, 1; -1i, 1i, -1, 1; 0, 1, 0, 2; 0, 1i, 0, 2i]);
        initVec = ones(4,1) + 1i * ones(4,1);
    end

    methods (Test)
        function testThrowBadHermiticity(tc)
            tc.fatalAssertError(@()EigenSolver(tc.effL, true, ...
                        tc.initVec, 1E-9), 'EigenSolver:badHermiticity');
        end

        function testHermitian(tc)
            HERMITIAN = true;
            H = ctranspose(tc.effL) * tc.effL;
            [vec, val] = EigenSolver(H, HERMITIAN, tc.initVec);

            actualVal = 0;

            actualVec = [-0.450257651351328 - 0.088285690737856i;
                         -0.719529322395919 + 0.092700289568880i;
                          0.178778680891222 - 0.315621815829032i;
                          0.359764661197959 - 0.046350144784440i];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - vec), tc.absTol);
        end

        function testNonHermitian(tc)
            HERMITIAN = false;
            [vec, val] = EigenSolver(tc.effL, HERMITIAN, tc.initVec);

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
