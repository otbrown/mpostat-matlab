% StationaryTest.m
% note that these tests should test the control flow of Stationary, but not
% whether or not it finds the correct answers
% Oliver Thomson Brown
% 2016-05-09

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true), matlab.unittest.fixtures.PathFixture('../../external/primme/Matlab')}) StationaryTest < matlab.unittest.TestCase

    properties
        absTol = 1E-11;
        sampleSz = 200;
        dmpoStat;
        dmpoStatH;
        eigTrack;
        eigTrackH;
        HILBY = 3;
        LENGTH = 5;
        COMPRESS = 27;
        dmpoInit;
        mpo;
        THRESHOLD = 1E-11;
    end

    methods (TestClassSetup)
        function ClassSetup(tc)
            % initialise dmpo
            tc.dmpoInit = MixDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);

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
                                     [tc.HILBY, tc.HILBY, tc.HILBY, ...
                                     tc.HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, ...
                                     tc.HILBY]);
            lmpo(:, :, :, :, 2, 2) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, ...
                                     tc.HILBY]);

            tc.mpo = cell(tc.LENGTH, 1);
            tc.mpo{1} = lmpo(:, :, :, :, 2, :);
            tc.mpo{tc.LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                tc.mpo{site} = lmpo;
            end
            mpoH = MPOHermProd(tc.mpo);

            % solve using Stationary
            [tc.dmpoStat, tc.eigTrack] = Stationary(tc.dmpoInit, ...
                                         tc.mpo, tc.THRESHOLD, 'direct');
            [tc.dmpoStatH, tc.eigTrackH] = Stationary(tc.dmpoInit,  ...
                                           mpoH, tc.THRESHOLD, ...
                                           'hermitian');
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.dmpoStat, 'cell');
            tc.fatalAssertClass(tc.dmpoStatH, 'cell');
            tc.fatalAssertClass(tc.eigTrack, 'double');
            tc.fatalAssertClass(tc.eigTrackH, 'double');
        end

        function testThrowBadHERMITIAN(tc)
            BAD_HERMITIAN = 'test';
            tc.fatalAssertError(@()Stationary(tc.dmpoInit, tc.mpo, ...
                                    tc.THRESHOLD, BAD_HERMITIAN), ...
                                    'Stationary:badHERMITIAN');
        end

        function testThrowBadHermiticity(tc)
            tc.fatalAssertError(@()Stationary(tc.dmpoInit, tc.mpo, ...
                                    tc.THRESHOLD, 'hermitian'), ...
                                    'EigenSolver:badHermiticity');
        end

        function testShape(tc)
            tc.fatalAssertSize(tc.dmpoStat, [tc.LENGTH, 1]);
            tc.fatalAssertSize(tc.dmpoStatH, [tc.LENGTH, 1]);
            tc.fatalAssertSize(tc.eigTrack, [2*(tc.LENGTH-1), 1]);
            tc.fatalAssertSize(tc.eigTrackH, [2*(tc.LENGTH-1), 1]);
        end

        function testTrace(tc)
            tr = DMPOTrace(tc.dmpoStat);
            trH = DMPOTrace(tc.dmpoStatH);
            tc.assertLessThan((tr-1), tc.absTol);
            tc.assertLessThan((trH-1), tc.absTol);
        end

        function testEigZero(tc)
            tc.assertLessThan(abs(tc.eigTrack(end)), tc.THRESHOLD);
            tc.assertLessThan(abs(tc.eigTrackH(end)), tc.THRESHOLD);
        end

        function testZZZ(tc)
            % since we are using a purely dissipative Liouvillian, we
            % expect the final state to be |000><000|, so we check that
            % element is 1
            zzz = tc.dmpoStat{1}(:, :, 1, 1)  ...
                    * tc.dmpoStat{2}(:, :, 1, 1) ...
                    * tc.dmpoStat{3}(:, :, 1, 1);
            zzzH = tc.dmpoStatH{1}(:, :, 1, 1) ...
                    * tc.dmpoStatH{2}(:, :, 1, 1) ...
                    * tc.dmpoStatH{3}(:, :, 1, 1);
            tc.assertLessThan((zzz - 1), tc.absTol);
            tc.assertLessThan((zzzH - 1), tc.absTol);
        end

        function testZeroes(tc)
            % since the density matrix should be zero everywhere but the
            % first element, we test a sample of them
            SPACE = tc.HILBY^tc.LENGTH;
            for testNum = 1 : 1 : tc.sampleSz
                % select random number starting at 1 for braState to avoid
                % getting the all zeroes state
                braState = randi([1, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH);
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH);

                coefft = 1;
                coefftH = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site) + 1;
                    ket = braBits(site) + 1;
                    coefft = coefft * tc.dmpoStat{site}(:, :, bra, ket);
                    coefftH = coefftH * tc.dmpoStatH{site}(:, :, bra, ket);
                end
                tc.assertLessThan(coefft, tc.absTol);
                tc.assertLessThan(coefftH, tc.absTol);
            end
        end
    end
end
