% StationaryTest.m
% note that these tests should test the control flow of Statioanry, but not
% whether or not it finds the correct answers
% Oliver Thomson Brown
% 2016-05-09

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) StationaryTest < matlab.unittest.TestCase

    properties
        dmpoStat;
        eigTrack;
        HILBY = 3;
        LENGTH = 3;
        COMPRESS = 0;
        dmpoInit;
        mpo;
        THRESHOLD = 1E-7;
        RUNMAX = 100;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            % initialise dmpo
            tc.dmpoInit = SuperDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);

            % build liouvillian mpo for dissipation only
            gamma = 1;
            I = eye(tc.HILBY);
            exDwn = zeros(tc.HILBY);
            for i = 1 : 1 : (tc.HILBY - 1)
                exDwn(i, i+1) = sqrt(i);
            end

            ident = kron(I, I);
            dissipator = 0.5 * gamma * (kron(2*conj(exDwn), exDwn) ...
                         - kron(I, ctranspose(exDwn)*exDwn) ...
                         - kron(transpose(exDwn)*conj(exDwn), I));

            lmpo = zeros(tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, 2, 2);
            lmpo(:, :, :, :, 2, 1) = reshape(dissipator, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 2, 2) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            tc.mpo = cell(tc.LENGTH, 1);
            tc.mpo{1} = lmpo(:, :, :, :, 2, :);
            tc.mpo{tc.LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                tc.mpo{site} = lmpo;
            end

            % solve using Stationary
            [tc.dmpoStat, tc.eigTrack] = Stationary(tc.dmpoInit, tc.mpo, ...
                                                    tc.THRESHOLD, tc.RUNMAX);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.dmpoStat, 'cell');
            tc.fatalAssertClass(tc.eigTrack, 'double');
        end

        function testShape(tc)
            tc.fatalAssertSize(tc.dmpoStat, [tc.LENGTH, 1]);
            tc.fatalAssertSize(tc.eigTrack, [tc.RUNMAX-1, 1]);
        end
    end
end
