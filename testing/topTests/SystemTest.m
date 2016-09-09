% SystemTest.m
% tests the variational stationary state search on some systems with known
% results
% Oliver Thomson Brown
% 2016-05-12

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) SystemTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        dispTol = 1E-4;
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
            % no Hamiltonian, just a dissipative turn -- trivial
            % stationary state
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
            [dmpoStat, eigTrack] = Stationary(tc.dmpoInit, mpo, ...
                                              tc.THRESHOLD, tc.RUNMAX);

			eigTrack = eigTrack(~isnan(eigTrack));

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
                if braState == 0 && ketState == 0
                    braState = 1;
                    ketState = 1;
                end
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
            tc.assertLessThan(abs(tr - 1), tc.absTol);
            tc.assertLessThan(abs(zzz - 1), tc.absTol);
            tc.assertLessThan(abs(eigTrack(end)), tc.THRESHOLD);
            for sample = 1 : 1 : sampleSize
                tc.assertLessThan(abs(elements{sample}), tc.absTol);
            end
        end

        function nonDriven(tc)
            % the Hamiltonian from the commensurate density investigation
            % minus the crucial driving term -- still has a trivial stationary
            % state
            % build Liouvillian mpo
            % parameters
            detuning01 = 0;
            detuning02 = 0;
            hop = 0.5;
            diss21 = 1;
            diss10 = 0.1;

            % local operators
            I = eye(tc.HILBY);
            a = zeros(tc.HILBY);
            for i = 1 : 1 : (tc.HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end
            a10 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            a21 = [0, 0, 0; 0, 0, sqrt(2); 0, 0, 0];

            % mpo building blocks
            ident = kron(I, I);
            Hloc = [0, 0, 0; 0, detuning01, 0; 0, 0, detuning02];
            dissipator = 0.5 * diss10 * (kron(2 * conj(a10), a10) ...
                         - kron(I, ctranspose(a10) * a10) ...
                         - kron(transpose(a10) * conj(a10), I)) ...
                         + ...
                         0.5 * diss21 * (kron(2 * conj(a21), a21) ...
                         - kron(I, ctranspose(a21) * a21) ...
                         - kron(transpose(a21) * conj(a21), I));

            % put liouvillian mpo together
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo = zeros(tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, 6, 6);
            lmpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 6) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*hop*ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 3) = reshape(kron(-1i*hop*ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*hop*a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 5) = reshape(kron(-1i*hop*a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            % build mpo cell
            mpo = cell(tc.LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 6, :);
            mpo{tc.LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(tc.dmpoInit, mpo, ...
                                              tc.THRESHOLD, tc.RUNMAX);
			eigTrack = eigTrack(~isnan(eigTrack));

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
                if braState == 0 && ketState == 0
                    braState = 1;
                    ketState = 1;
                end
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
            tc.assertLessThan(abs(tr - 1), tc.absTol);
            tc.assertLessThan(abs(zzz - 1), tc.absTol);
            tc.assertLessThan(abs(eigTrack(end)), tc.THRESHOLD);
            for sample = 1 : 1 : sampleSize
                tc.assertLessThan(abs(elements{sample}), tc.absTol);
            end
        end

        function drivenDiss3Site(tc)
            % the Hamiltonian from the commensurate density investigation
            % build Liouvillian mpo
            % parameters
            detuning01 = 0;
            detuning02 = 0;
            hop = 0.5;
            las02Intensity = 8;
            diss21 = 1;
            diss10 = 0.1;

            % local operators
            I = eye(tc.HILBY);
            a = zeros(tc.HILBY);
            for i = 1 : 1 : (tc.HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end
            a10 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            a21 = [0, 0, 0; 0, 0, sqrt(2); 0, 0, 0];

            % mpo building blocks
            ident = kron(I, I);
            Hloc = [0, 0, conj(las02Intensity); ...
                    0, detuning01, 0; ...
                    las02Intensity, 0, detuning02];
            dissipator = 0.5 * diss10 * (kron(2 * conj(a10), a10) ...
                         - kron(I, ctranspose(a10) * a10) ...
                         - kron(transpose(a10) * conj(a10), I)) ...
                         + ...
                         0.5 * diss21 * (kron(2 * conj(a21), a21) ...
                         - kron(I, ctranspose(a21) * a21) ...
                         - kron(transpose(a21) * conj(a21), I));

            % put liouvillian mpo together
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo = zeros(tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, 6, 6);
            lmpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 6) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*hop*ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 3) = reshape(kron(-1i*hop*ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*hop*a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:, :, :, :, 6, 5) = reshape(kron(-1i*hop*a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            % build mpo cell
            mpo = cell(tc.LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 6, :);
            mpo{tc.LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(tc.dmpoInit, mpo, ...
                                              tc.THRESHOLD, tc.RUNMAX);
			eigTrack = eigTrack(~isnan(eigTrack));

            % calculate some values to assert against
            tr = DMPOTrace(dmpoStat);
            rho000 = 1;
            rho111 = 1;
            rho222 = 1;
            for site = 1 : 1 : tc.LENGTH
                rho000 = rho000 * dmpoStat{site}(:, :, 1, 1);
                rho111 = rho111 * dmpoStat{site}(:, :, 2, 2);
                rho222 = rho222 * dmpoStat{site}(:, :, 3, 3);
            end

            % make some assertions
            tc.assertLessThan(abs(tr - 1), tc.absTol);
            tc.assertLessThan(abs(eigTrack(end)), tc.THRESHOLD);
            tc.assertLessThan(abs(rho000 - 5.665E-4), tc.dispTol);
            tc.assertLessThan(abs(rho111 - 0.7377), tc.dispTol);
            tc.assertLessThan(abs(rho222 - 1.3439E-4), tc.dispTol);
        end
    end

end
