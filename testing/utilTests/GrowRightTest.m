% GrowRightTest.m
% Oliver Thomson Brown
% 2016-03-22

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) GrowRightTest < matlab.unittest.TestCase

    properties
        absTol = 1E-15;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        dmpo;
        rhoVec;
        iRight;
        lRight;
        impo;
        lmpo;
        L;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 3};
        testLENGTH = {3, 4, 3};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.dmpo = DDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.dmpo = Can(tc.dmpo, [tc.LENGTH : -1 : 2], 'R');
            tc.rhoVec = ones(tc.HILBY^(2*tc.LENGTH), 1) / (tc.HILBY^tc.LENGTH);
            tc.iRight = cell(tc.LENGTH, 1);
            tc.lRight = cell(tc.LENGTH, 1);
            tc.iRight{tc.LENGTH} = 1;
            tc.lRight{tc.LENGTH} = 1;
            tc.impo = reshape(eye(tc.HILBY^2), [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            tc.lmpo = cell(tc.LENGTH, 1);

            % build liouvillian mpo
            mpo = zeros(tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, 6, 6);

            U = 1;
            J = 0.5;
            gamma = 0.1;

            I = eye(tc.HILBY);
            a = zeros(tc.HILBY);
            for i = 1 : 1 : (tc.HILBY - 1)
                a(i, i + 1) = sqrt(i);
            end

            ident = kron(I, I);
            Hloc = zeros(tc.HILBY);
            for i = 2 : 1 : tc.HILBY
                Hloc(i, i) = i * U;
            end
            dissipator = 0.5 * gamma * (kron(2 * conj(a), a) ...
                         - kron(I, ctranspose(a) * a) ...
                         - kron(transpose(a) * conj(a), I));
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            mpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 1, 1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 6, 6) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*J*ctranspose(a)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 6, 3) = reshape(kron(-1i*J*ctranspose(a), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*J*a), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            mpo(:, :, :, :, 6, 5) = reshape(kron(-1i*J*a, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            tc.lmpo{1} = mpo(:, :, :, :, 6, :);
            tc.lmpo{tc.LENGTH} = mpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                tc.lmpo{site} = mpo;
            end

            % populate rightBlock
            for site = tc.LENGTH : -1 : 2
                [rowSz, colSz, ~, ~] = size(tc.dmpo{site});
                [~, ~, ~, ~, OP_ROW, OP_COL] = size(tc.lmpo{site});
                tc.iRight{site - 1} = GrowRight(tc.dmpo{site}, tc.impo, tc.iRight{site}, rowSz, colSz, tc.HILBY, 1);
                tc.lRight{site - 1} = GrowRight(tc.dmpo{site}, tc.lmpo{site}, tc.lRight{site}, rowSz, colSz, tc.HILBY, OP_ROW);
            end

            % build Liouvillian
            IH = eye(tc.HILBY^tc.LENGTH);
            a_n = zeros(tc.HILBY^tc.LENGTH, tc.HILBY^tc.LENGTH, tc.LENGTH);
            H_n = zeros(tc.HILBY^tc.LENGTH, tc.HILBY^tc.LENGTH, tc.LENGTH);
            for target = 1 : 1 : tc.LENGTH
                op = 1;
                ham = 1;
                for site = 1 : 1 : (target - 1)
                    op = kron(op, I);
                    ham = kron(ham, I);
                end
                op = kron(op, a);
                ham = kron(ham, Hloc);
                for site = (target + 1) : 1 : tc.LENGTH
                    op = kron(op, I);
                    ham = kron(ham, I);
                end
                a_n(:,:,target) = op;
                H_n(:,:,target) = ham;
            end

            Hint = zeros(tc.HILBY^tc.LENGTH);
            for site = 1 : 1 : (tc.LENGTH - 1)
                Hint = Hint + a_n(:,:,site) * ctranspose(a_n(:,:,site + 1)) ...
                            + ctranspose(a_n(:,:,site)) * a_n(:,:,site+1);
            end
            Hint = - J * Hint;

            H = zeros(tc.HILBY^tc.LENGTH);
            for site = 1 : 1 : tc.LENGTH
                H = H + H_n(:, :, site);
            end
            H = H + Hint;

            tc.L = -1i * (kron(IH, H) - kron(H, IH));
            for site = 1 : 1 : tc.LENGTH
                tc.L = tc.L + 0.5 * gamma ...
                * (kron(2*conj(a_n(:,:,site)), a_n(:,:,site)) ...
                - kron(IH, ctranspose(a_n(:,:,site))*a_n(:,:,site)) ...
                - kron(transpose(a_n(:,:,site))*a_n(:,:,site), IH));
            end
        end
    end

    methods (Test)
        function vecNormTest(tc)
            % calculate vector norm by using LeftGrow past the end of the system
            [rowSz, colSz, ~, ~] = size(tc.dmpo{1});
            testVecNorm = GrowRight(tc.dmpo{1}, tc.impo, tc.iRight{1}, rowSz, colSz, tc.HILBY, 1);

            % calculate the vector norm the old fashioned way...
            vecNorm = ctranspose(tc.rhoVec) * tc.rhoVec;

            epsilon = abs(testVecNorm - vecNorm);
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function liouvExpTest(tc)
            % calculate Liouvillian expectation by growin g past the last site
            [rowSz, colSz, ~, ~] = size(tc.dmpo{1});
            testRExp = GrowRight(tc.dmpo{1}, tc.lmpo{1}, tc.lRight{1}, rowSz, colSz, tc.HILBY, 1);

            % calculate the Liouvillian expectation the old-fashioned way
            lExp = ctranspose(tc.rhoVec) * tc.L * tc.rhoVec;

            epsilon = abs(testRExp - lExp);
            tc.assertLessThan(epsilon, tc.absTol);
        end
    end
end
