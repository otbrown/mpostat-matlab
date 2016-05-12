% SystemTest.m
% tests the variational stationary state search on some systems with known
% results
% Oliver Thomson Brown
% 2016-05-12

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) SystemTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        HILBY = 3;
        LENGTH = 3;
        COMPRESS = 100;
        dmpoInit;
        THRESHOLD = 1E-7;
        RUNMAX = 50;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            tc.dmpoInit = SuperDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
        end
    end

    methods (Test)
        function dissipatorOnly(tc)
            % build Liouvillian mpo
            gamma = 1;
            I = eye(tc.HILBY);
            a = zeros(tc.HILBY);
            for i = 1 : 1 : (tc.HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end

            ident = kron(I, I);
            dissipator = 0.5 * gamma * (kron(2 * conj(a), a) ...
                         - kron(I, ctranspose(a) * a) ...
                         - kron(transpose(a) * conj(a), I));

            lmpo = zeros(tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, 2, 2);
            lmpo(:, :, :, :, 2, 1) = reshape(dissipator, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 2, 2) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            mpo = cell(tc.LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 2, :);
            mpo{tc.LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(tc.dmpoInit, mpo, tc.COMPRESS, ...
                                              tc.THRESHOLD, tc.RUNMAX);

            % calculate some values to assert against
            SPACE = tc.HILBY ^ tc.LENGTH;
            sampleSize = min(floor(0.25 * SPACE^2), 100);
            elements = cell(sampleSize, 1);
            [elements{:}] = deal(1);

            tr = DMPOTrace(dmpoStat);

            zzz = 1;
            for site = 1 : 1 : tc.LENGTH
                zzz = zzz * dmpoStat{site}(:, :, 1, 1);
            end

            for sample = 1 : 1 : sampleSize
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH);
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH);

                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site) + 1;
                    ket = ketBits(site) + 1;
                    elements{sample} = elements{sample} ...
                                       * dmpoStat{site}(:, :, bra, ket);
                end
            end

            % make some assertions
            tc.assertLessThan((tr - 1), tc.absTol);
            tc.assertLessThan((zzz - 1), tc.absTol);
            tc.assertLessThan(eigTrack(end), tc.absTol);
            for sample = 1 : 1 : sampleSize
                tc.assertLessThan(elements{sample}, tc.absTol);
            end
        end
    end

end
