% GrowBlockTest.m
% Oliver Thomson Brown
% 2016-10-24

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) GrowBlockTest < matlab.unittest.TestCase

    properties
        HILBY = 2;
        LENGTH = 4;
        COMPRESS = 0;
        dmpo;
        mpo;
        initLeft;
        initRight;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.initLeft = cell(tc.LENGTH, 1);
            tc.initRight = cell(tc.LENGTH, 1);

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

            tc.mpo{1} = mpo(:, :, :, :, 6, :);
            tc.mpo{tc.LENGTH} = mpo(:, :, :, :, :, 1);
            for site = 2 : 1 : (tc.LENGTH - 1)
                tc.mpo{site} = mpo;
            end
        end
    end

    methods (Test)
        function testThrowBadDirection(tc)
            badDirection = 'm';
            tc.fatalAssertError(@()GrowBlock(tc.dmpo, tc.mpo, tc.initLeft, ...
                                            tc.initRight, 1, badDirection), ...
                                            'GrowBlock:BadDirection');
        end

        function testGrowLeft(tc)
            direction = 'L';
            left = tc.initLeft;
            left{1} = 1;
            for site = 2 : 1 : tc.LENGTH
                left{site} = GrowBlock(tc.dmpo, tc.mpo, left, tc.initRight, ...
                                        (site-1), direction);
            end

            % be assertive!
            for site = 1 : 1 : tc.LENGTH
                tc.assertClass(left{site}, 'double');
                tc.assertNotEmpty(left{site});
            end
        end

        function testGrowRight(tc)
            direction = 'R';
            right = tc.initRight;
            right{tc.LENGTH} = 1;
            for site = (tc.LENGTH - 1) : -1 : 1
                right{site} = GrowBlock(tc.dmpo, tc.mpo, tc.initLeft, ...
                                        right, (site+1), direction);
            end

            % be assertive!
            for site = 1 : 1 : tc.LENGTH
                tc.assertClass(right{site}, 'double');
                tc.assertNotEmpty(right{site});
            end
        end
    end
end
