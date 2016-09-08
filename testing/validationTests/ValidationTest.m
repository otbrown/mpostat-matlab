% ValidationTest.m
% validate mpostat code against known results
% Oliver Thomson Brown
% 2016-9-06

% This set of tests requires the following set of .mat files to be present:
%   - FourSiteExact.mat
%   - ElevenSiteTEBD.mat
%   - MFS15.mat
% These files contain the results we validate against.

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) ValidationTest < matlab.unittest.TestCase

    properties
        THRESHOLD = 1E-7;
        RUNMAX = 500;
        COMPRESS = 100;
    end

    methods (Test)
        function FourSiteExact(tc)
            % system parameters
            HILBY = 3;
            LENGTH = 4;
            hop = 0.5;
            las02Intensity = 8;
            detuning01 = 0;
            detuning02 = 0;
            diss21 = 1;
            diss10 = 0.1;
            fname = 'data/FourSiteExact.mat';

            % inputs
            dmpoInit = SuperDMPO(HILBY, LENGTH, tc.COMPRESS);

            % local operators
            I = eye(HILBY);
            a = zeros(HILBY);
            for i = 1 : 1 : (HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end
            a10 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            a21 = [0, 0, 0; 0, 0, sqrt(2); 0, 0, 0];

            % mpo building blocks
            ident = kron(I, I);
            Hloc = [0, 0, conj(las02Intensity) / sqrt(2); ...
                    0, detuning01, 0; ...
                    las02Intensity / sqrt(2), 0, detuning02];
            dissipator = 0.5 * diss10 * (kron(2 * conj(a10), a10) ...
                         - kron(I, ctranspose(a10) * a10) ...
                         - kron(transpose(a10) * conj(a10), I)) ...
                         + ...
                         0.5 * diss21 * (kron(2 * conj(a21), a21) ...
                         - kron(I, ctranspose(a21) * a21) ...
                         - kron(transpose(a21) * conj(a21), I));

            % put liouvillian mpo together
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo = zeros(HILBY, HILBY, HILBY, HILBY, 6, 6);
            lmpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 6) = reshape(ident, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*hop*ctranspose(a)), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 3) = reshape(kron(-1i*hop*ctranspose(a), I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*hop*a), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 5) = reshape(kron(-1i*hop*a, I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);

            % build mpo cell
            mpo = cell(LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 6, :);
            mpo{LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, tc.THRESHOLD, tc.RUNMAX);
            eigTrack = eigTrack(~isnan(eigTrack));

            % calculate observables
            % allocate
            testN = zeros(LENGTH,1);
            testAdAdAA = zeros(LENGTH,1);
            testG1 = zeros(LENGTH);
            testG2 = zeros(LENGTH);

            idOp = zeros(HILBY, HILBY, LENGTH);
            for site = 1 : 1 : LENGTH
                idOp(:, :, site) = eye(HILBY);
            end
            op = idOp;

            for site = 1 : 1 : LENGTH
                op(:, :, site) = ctranspose(a) * a;
                testN(site) = DMPOExp(dmpoStat, op);
                op = idOp;
                op(:, :, site) = ctranspose(a)*ctranspose(a)*a*a;
                testAdAdAA(site) = DMPOExp(dmpoStat, op);
                op = idOp;
            end

            for j = 1 : 1 : LENGTH
                for l = 1 : 1 : LENGTH
                    if j ~= l
                        op(:, :, j) = ctranspose(a);
                        op(:, :, l) = a;
                        expAdA = DMPOExp(dmpoStat, op);
                        op = idOp;
                        op(:, :, j) = ctranspose(a) * a;
                        op(:, :, l) = ctranspose(a) * a;
                        expAdAdAA = DMPOExp(dmpoStat, op);
                        op = idOp;
                    else
                        expAdA = Nmpo(j);
                        op(:, :, j) = ctranspose(a)*ctranspose(a)*a*a;
                        expAdAdAA = DMPOExp(dmpoStat, op);
                        op = idOp;
                    end
                    testG1(j,l) = expAdA / sqrt(Nmpo(j)*Nmpo(l));
                    testG2(j,l) = expAdAdAA / (Nmpo(j)*Nmpo(l));
                end
            end

            % be assertive!
            load(fname);
            tc.assertLessThan(abs(N - testN), tc.THRESHOLD);
            tc.assertLessThan(abs(AdAdAA - testAdAdAA), tc.THRESHOLD);
            tc.assertLessThan(abs(G1 - testG1), tc.THRESHOLD);
            tc.assertLessThan(abs(G2 - testG2), tc.THRESHOLD);
        end

        function ElevenSiteTEBD(tc)
            % system parameters
            HILBY = 3;
            LENGTH = 11;
            hop = 0.5;
            las02Intensity = 8;
            detuning01 = 0;
            detuning02 = 0;
            diss21 = 1;
            diss10 = 0.1;
            fname = 'data/ElevenSiteTEBD.mat';

            % inputs
            dmpoInit = SuperDMPO(HILBY, LENGTH, tc.COMPRESS);

            % local operators
            I = eye(HILBY);
            a = zeros(HILBY);
            for i = 1 : 1 : (HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end
            a10 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            a21 = [0, 0, 0; 0, 0, sqrt(2); 0, 0, 0];

            % mpo building blocks
            ident = kron(I, I);
            Hloc = [0, 0, conj(las02Intensity) / sqrt(2); ...
                    0, detuning01, 0; ...
                    las02Intensity / sqrt(2), 0, detuning02];
            dissipator = 0.5 * diss10 * (kron(2 * conj(a10), a10) ...
                         - kron(I, ctranspose(a10) * a10) ...
                         - kron(transpose(a10) * conj(a10), I)) ...
                         + ...
                         0.5 * diss21 * (kron(2 * conj(a21), a21) ...
                         - kron(I, ctranspose(a21) * a21) ...
                         - kron(transpose(a21) * conj(a21), I));

            % put liouvillian mpo together
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo = zeros(HILBY, HILBY, HILBY, HILBY, 6, 6);
            lmpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 6) = reshape(ident, ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*hop*ctranspose(a)), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 3) = reshape(kron(-1i*hop*ctranspose(a), I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*hop*a), ...
                                     [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 6, 5) = reshape(kron(-1i*hop*a, I), ...
                                     [HILBY, HILBY, HILBY, HILBY]);

            % build mpo cell
            mpo = cell(LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 6, :);
            mpo{LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, tc.THRESHOLD, tc.RUNMAX);
            eigTrack = eigTrack(~isnan(eigTrack));

            % calculate observables
            % allocate
            testSite = 5;
            testG1 = zeros(LENGTH);
            testG2 = zeros(LENGTH);

            idOp = zeros(HILBY, HILBY, LENGTH);
            for site = 1 : 1 : LENGTH
                idOp(:, :, site) = eye(HILBY);
            end
            op = idOp;

            op(:, :, testSite) = ctranspose(a)*a;
            testN = DMPOExp(dmpoStat, op);
            op = idOp;

            op(:, :, testSite) = ctranspose(a)*ctranspose(a)*a*a;
            testAdAdAA = DMPOExp(dmpoStat, op);
            op = idOp;

            for j = 1 : 1 : LENGTH
                for l = 1 : 1 : LENGTH
                    if j ~= l
                        op(:, :, j) = ctranspose(a);
                        op(:, :, l) = a;
                        expAdA = DMPOExp(dmpoStat, op);
                        op = idOp;
                        op(:, :, j) = ctranspose(a) * a;
                        op(:, :, l) = ctranspose(a) * a;
                        expAdAdAA = DMPOExp(dmpoStat, op);
                        op = idOp;
                    else
                        expAdA = Nmpo(j);
                        op(:, :, j) = ctranspose(a)*ctranspose(a)*a*a;
                        expAdAdAA = DMPOExp(dmpoStat, op);
                        op = idOp;
                    end
                    testG1(j,l) = expAdA / sqrt(Nmpo(j)*Nmpo(l));
                    testG2(j,l) = expAdAdAA / (Nmpo(j)*Nmpo(l));
                end
            end

            % be assertive!
            load(fname);
            tc.assertLessThan(abs(N - testN), tc.THRESHOLD);
            tc.assertLessThan(abs(AdAdAA - testAdAdAA), tc.THRESHOLD);
            tc.assertLessThan(abs(G1 - testG1), tc.THRESHOLD);
            tc.assertLessThan(abs(G2 - testG2), tc.THRESHOLD);
        end

        function MFS15(tc)
            % system parameters
            HILBY = 2;
            LENGTH = 15;
            J = 1;
            h = J;
            V = 0.5 * J;
            diss = 0.1 * J;
            fname = 'data/MFS15.mat';

            % inputs
            dmpoInit = SuperDMPO(HILBY, LENGTH, tc.COMPRESS);

            % local operators
            I = eye(HILBY);
            X = [0, 1; 1, 0];
            Y = [0, -1i; 1i, 0];
            Z = [1, 0; 0, -1];

            % mpo building blocks
            ident = kron(I, I);
            Hloc = h * Z;
            dissipator = kron(2*conj(K), K) ...
             - kron(I, ctranspose(K)*K) - kron(transpose(K)*conj(K), I);

            % put liouvillian mpo together
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                            [HILBY, HILBY,HILBY, HILBY]);
            lmpo(:, :, :, :, 2, 1) = reshape(kron(I, X), ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 3, 1) = reshape(kron(X, I), ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 1) = reshape(mpoLoc,
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 2) = reshape(kron(I, -1i*J*X), ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 3) = reshape(kron(1i*J*X, I), ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 4, 4) = reshape(ident, ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 2, 2) = reshape(V * ident / J, ...
                                            [HILBY, HILBY, HILBY, HILBY]);
            lmpo(:, :, :, :, 3, 3) = reshape(V * ident / J, ...
                                            [HILBY, HILBY, HILBY, HILBY]);

            % build mpo cell
            mpo = cell(LENGTH, 1);
            mpo{1} = lmpo(:, :, :, :, 4, :);
            mpo{LENGTH} = lmpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (LENGTH - 1)
                mpo{site} = lmpo;
            end

            % solve using Stationary
            [dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, tc.THRESHOLD, tc.RUNMAX);
            eigTrack = eigTrack(~isnan(eigTrack));

            % calculate observables
            testSite = 8;
            testCorXX = zeros(4, 1);
            XOp = zeros(HILBY, HILBY, LENGTH);
            for site = 1 : 1 : LENGTH
                XOp(:, :, site) = I;
            end
            XOp(:, :, testSite) = X;
            op = XOp;

            for l = 1 : 1 : 4;
                op(:, :, testSite + l) = X;
                testCorXX(l) = DMPOExp(dmpoStat, op);
                op = XOp;
            end

            % be assertive!
            load(fname);
            for l = 1 : 1 : 4
                tc.assertLessThan(abs(CorXX(l) - testCorXX(l)), tc.THRESHOLD);
            end
        end
end