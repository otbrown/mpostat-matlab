% MPOHermProdTest.m
% Oliver Thomson Brown
% 2016-10-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) MPOHermProdTest < matlab.unittest.TestCase

    properties
        HILBY = 2;
        LENGTH = 4;
        OPDIM = 6;
        absTol = 1E-14;
        mpo;
        hpMPO;
        liouv;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            % build local operators
            DIM = tc.HILBY^tc.LENGTH;
            I = eye(tc.HILBY);
            Is = eye(DIM);

            A = zeros(tc.HILBY);
            for index = 1 : 1 : tc.HILBY-1
                A(index, index+1) = sqrt(index) * (1 + 1i)/sqrt(2);
            end

            Hloc = ctranspose(A)*A + ctranspose(A) + A;

            % build Hamiltonian
            siteOp = zeros(DIM, DIM, tc.LENGTH);
            for target = 1 : 1 : tc.LENGTH
                op = 1;
                for site = 1 : 1 : (target - 1)
                    op = kron(op, I);
                end
                op = kron(op, A);
                for site = (target + 1) : 1 : tc.LENGTH
                    op = kron(op, I);
                end
                siteOp(:, :, target) = op;
            end

            local = zeros(DIM);
            for site = 1 : 1 : tc.LENGTH
                local = local ...
                        + ctranspose(siteOp(:,:,site))*siteOp(:,:,site) ...
                        + ctranspose(siteOp(:,:,site)) + siteOp(:,:,site);
            end

            interaction = zeros(DIM);
            for site = 1 : 1 : (tc.LENGTH - 1)
                fwd = siteOp(:,:,site)*ctranspose(siteOp(:,:,site+1));
                rev = ctranspose(siteOp(:,:,site))*siteOp(:,:,site+1);
                interaction = interaction + fwd + rev;
            end

            H = local + interaction;

            % build Liouvillian
            siteDiss = zeros(DIM^2, DIM^2, tc.LENGTH);
            for site = 1 : 1 :tc.LENGTH
                op = siteOp(:, :, site);
                siteDiss(:, :, site) = 0.5 * (kron(2*conj(op), op) ...
                                        - kron(Is, ctranspose(op)*op) ...
                                        - kron(transpose(op)*conj(op), Is));
            end

            L = -1i * (kron(Is, H) - kron(H, Is));
            for site = 1 : 1 : tc.LENGTH
                L = L + siteDiss(:, :, site);
            end
            tc.liouv = L;

            % build mpo
            ident = kron(I, I);
            dissipator = 0.5 * (kron(2*conj(A), A) ...
                            - kron(I, ctranspose(A)*A) ...
                            - kron(transpose(A)*conj(A), I));
            mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

            lmpo = zeros(tc.HILBY,tc.HILBY,tc.HILBY,tc.HILBY,tc.OPDIM,tc.OPDIM);
            lmpo(:,:,:,:,tc.OPDIM,1) = reshape(mpoLoc, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,1,1) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,tc.OPDIM,tc.OPDIM) = reshape(ident, ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,tc.OPDIM,2) = reshape(kron(I, -1i*A), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,tc.OPDIM,3) = reshape(kron(1i*A, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,tc.OPDIM,4) = reshape(kron(I, -1i*ctranspose(A)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,tc.OPDIM,5) = reshape(kron(1i*ctranspose(A), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,2,1) = reshape(kron(I, ctranspose(A)), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,3,1) = reshape(kron(ctranspose(A), I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,4,1) = reshape(kron(I, A), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);
            lmpo(:,:,:,:,5,1) = reshape(kron(A, I), ...
                                     [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY]);

            tc.mpo = cell(tc.LENGTH, 1);
            tc.mpo{1} = lmpo(:,:,:,:,tc.OPDIM,:);
            tc.mpo{tc.LENGTH} = lmpo(:,:,:,:,:,1);
            for site = 2 : 1 : tc.LENGTH - 1
                tc.mpo{site} = lmpo;
            end

            % calculate hermitian product MPO
            tc.hpMPO = MPOHermProd(tc.mpo);
        end
    end

    methods (Test)
        function nonHermitianLiouv(tc)
            % check that the Liouvillian we've built is not Hermitian
            diff = abs(tc.liouv - ctranspose(tc.liouv));
            tc.fatalAssertTrue(any(diff(:)));
        end

        function systemSizeCheck(tc)
            tc.fatalAssertSize(tc.hpMPO, [tc.LENGTH, 1]);
        end

        function shapeCheck(tc)
            tc.fatalAssertSize(tc.hpMPO{1}, ...
                                [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, ...
                                1, tc.OPDIM^2]);
            for site = 2 : 1 : (tc.LENGTH - 1)
                tc.fatalAssertSize(tc.hpMPO{site}, ...
                                [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, ...
                                tc.OPDIM^2, tc.OPDIM^2]);
            end
            tc.fatalAssertSize(tc.hpMPO{tc.LENGTH}, ...
                                [tc.HILBY, tc.HILBY, tc.HILBY, tc.HILBY, ...
                                tc.OPDIM^2]);
        end

        function hermProdCheck(tc)
            SPACE = tc.HILBY^tc.LENGTH;

            % rebuild hermitian product of Liouvillian from mpo
            mpoHP = zeros(SPACE^2);
            for conjBra = 0 : 1 : (SPACE - 1)
                cBraBits = FWBase(conjBra, tc.HILBY, tc.LENGTH) + 1;
                for conjKet = 0 : 1 : (SPACE - 1)
                    cKetBits = FWBase(conjKet, tc.HILBY, tc.LENGTH) + 1;
                    lRow = conjKet * SPACE + conjBra + 1;
                    for bra = 0 : 1 : (SPACE - 1)
                        braBits = FWBase(bra, tc.HILBY, tc.LENGTH) + 1;
                        for ket = 0 : 1 : (SPACE - 1)
                            ketBits = FWBase(ket, tc.HILBY, tc.LENGTH) + 1;
                            lCol = ket * SPACE + bra + 1;

                            coefft = 1;
                            for site = 1 : 1 : tc.LENGTH
                                stateMPO = permute(tc.hpMPO{site}(cBraBits(site), cKetBits(site), braBits(site), ketBits(site), :, :), [5, 6, 1, 2, 3, 4]);
                                coefft = coefft * stateMPO;
                            end
                            mpoHP(lRow, lCol) = coefft;
                        end
                    end
                end
            end

            hP = ctranspose(tc.liouv) * tc.liouv;

            tc.assertLessThan(abs(ctranspose(mpoHP) - mpoHP), tc.absTol);
            tc.assertLessThan(abs(mpoHP - hP), tc.absTol);
        end
    end
end
