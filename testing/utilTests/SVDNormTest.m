% SVDNormTest.m
% Oliver Thomson Brown
% 2016-03-08

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) SVDNormTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        HILBY = 2;
        LENGTH = 7;
        COMPRESS = 100;
        dmpo;
        svdDMPO;
    end

    methods (TestMethodSetup)
        function MethodSetup(tc)
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.svdDMPO = SVDNorm(tc.dmpo);
        end
    end

    methods (Test)
        % Note that SVDNorm is in general an unsafe method on a DMPO, so it should
        % not be expected that it preserves trace, or any other property
        % virtual dimensions are also not preserved
        function testClass(tc)
            tc.fatalAssertClass(tc.svdDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.svdDMPO, size(tc.dmpo));
        end

        function testTensorShape(tc)
            % first site should still be 1 x C, physical dimensions should always
            % be HILBY
            colSz = size(tc.svdDMPO{1}, 2);
            tc.fatalAssertSize(tc.svdDMPO{1}, [1, colSz, tc.HILBY, tc.HILBY]);
            % last site should still be R x 1
            rowSz = size(tc.svdDMPO{tc.LENGTH}, 1);
            tc.fatalAssertSize(tc.svdDMPO{tc.LENGTH}, [rowSz, 1, tc.HILBY, tc.HILBY]);
            % can check physical dimensions on other sites
            for site = 2 : 1 : tc.LENGTH - 1
                [rowSz, colSz, ~, ~] = size(tc.svdDMPO{site});
                tc.fatalAssertSize(tc.svdDMPO{site}, [rowSz, colSz, tc.HILBY, tc.HILBY]);
            end
        end

        function leftNormTest(tc)
            % across all sites the new tensor is formed by reshaping the Q matrix
            % from a QR decomposition -- as a result A' * A = I (sort of)
            for site = 1 : 1 : tc.LENGTH
                [rowSz, colSz, ~, ~] = size(tc.svdDMPO{site});
                M = reshape(tc.svdDMPO{site}, [rowSz, colSz, tc.HILBY^2]);
                M = permute(M, [1, 3, 2]);
                M = reshape(M, [rowSz * tc.HILBY^2, colSz]);

                I = eye(colSz, colSz);
                epsilon = (ctranspose(M) * M) - I;
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end

        function vecNormTest(tc)
            % the SVD norm is used to normalise MPS which represent state vectors
            % so if we pretend dmpo is just a really big MPS we should find the
            % the rebuilt state vector rhoVec' * rhoVec = 1
            % the following code is "borrowed" from the debug function DMRebuild
            % obviously, be wary of using a large test system

            % RETURN ALLOCATION
            SPACE = tc.HILBY^tc.LENGTH;
            rhoVec = complex(zeros(SPACE^2, 1));

            % REBUILD VECTORISED DENSITY MATRIX
            for ketState = 0 : 1 : (SPACE - 1)
                for braState = 0 : 1 : (SPACE - 1)
                    stateDex = ketState * SPACE + braState + 1;
                    braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                    ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                    coefft = 1;
                    for site = 1 : 1 : tc.LENGTH
                        bra = braBits(site);
                        ket = ketBits(site);
                        coefft = coefft * tc.svdDMPO{site}(:, :, bra, ket);
                    end

                    rhoVec(stateDex) = coefft;
                end
            end

            vecNorm = ctranspose(rhoVec) * rhoVec;
            tc.assertLessThan((vecNorm - 1), tc.absTol);
        end
    end
end
