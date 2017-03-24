% EigenSolverTest.m
% Oliver Thomson Brown
% 2016-10-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true), matlab.unittest.fixtures.PathFixture('../../external/primme/Matlab')}) EigenSolverTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        effL = sparse([2, -1, -1, 1; -1i, 1i, -1, 1; 0, 1, 0, 2; 0, 1i, 0, 2i]);
        initVec = ones(4, 1);
    end

    methods (Test)
        function testThrowBadHermiticity(tc)
            tc.fatalAssertError(@()EigenSolver(tc.effL, true, false, tc.initVec, 1E-9), ...
            'EigenSolver:badHermiticity');
        end

        function testPrimme(tc)
            HERMITIAN = true;
            PRIMME = true;
            H = ctranspose(tc.effL) * tc.effL;
            [vec, val] = EigenSolver(H, HERMITIAN, PRIMME, tc.initVec, 1);

            actualVal = 0;
            actualVec = [ 0.302322027807821 + 0.345148819050046i;
                          0.626057451236754 + 0.366562214671157i;
                         -0.334442121239488 + 0.140454316093355i;
                         -0.313028725618377 - 0.183281107335578i];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - vec), tc.absTol);
        end

        function testHermitian(tc)
            HERMITIAN = true;
            PRIMME = false;
            H = ctranspose(tc.effL) * tc.effL;
            [vec, val] = EigenSolver(H, HERMITIAN, PRIMME, tc.initVec, 1);

            actualVal = 0;
            actualVec = [-0.269871683328030 + 0.371073564576041i;
                         -0.219270742704024 + 0.691546188528076i;
                         -0.210837252600023 - 0.295172153640033i;
                          0.109635371352012 - 0.345773094264038i];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - vec), tc.absTol);
        end

        function testNonHermitian(tc)
            HERMITIAN = false;
            PRIMME = false;
            [vec, val] = EigenSolver(tc.effL, HERMITIAN, PRIMME, ...
                                        tc.initVec);

            actualVal = 1.920620363439616 + 0.618180731783746i;
            actualVec = [0.909915536697958;
                         0.356041764971642;
                         0.150479197352232;
                         0.150479197352232];

            tc.assertLessThan(abs(actualVal - val), tc.absTol);
            tc.assertLessThan(abs(actualVec - abs(vec)), tc.absTol);
        end
    end
end
