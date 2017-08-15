% DMPOExpTest.m
% Oliver Thomson Brown
% 2016-03-11

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) DMPOExpTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        COMPRESS = 100;
        HILBY;
        LENGTH;
        STATE;
        dmpo;
        prodDMPO;
        idOp;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testLENGTH = {7, 6, 5, 4};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.STATE = tc.HILBY^tc.LENGTH - 1;
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.prodDMPO = ProdDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS, tc.STATE);
            tc.idOp = repmat(eye(tc.HILBY), [1, 1, tc.LENGTH]);
        end
    end

    methods (Test)
        function testClass(tc)
            expect = DMPOExp(tc.dmpo, tc.idOp);
            tc.fatalAssertClass(expect, 'double');
        end

        function testArbTrace(tc)
            idExp = DMPOExp(tc.dmpo, tc.idOp);
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(idExp) - abs(tr));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testProdTrace(tc)
            idExp = DMPOExp(tc.prodDMPO, tc.idOp);
            tr = DMPOTrace(tc.dmpo);
            epsilon = abs(abs(idExp) - abs(tr));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testNExpectation(tc)
            nSiteOp = cell(tc.LENGTH, 1);
            n = diag([0 : 1 : tc.HILBY - 1]);
            for site = 1 : 1 : tc.LENGTH
                nSiteOp{site} = repmat(eye(tc.HILBY), [1, 1, tc.LENGTH]);
                nSiteOp{site}(:, :, site) = n;
            end

            nExp = 0;
            for site = 1 : 1 : tc.LENGTH
                 nExp = nExp + DMPOExp(tc.prodDMPO, nSiteOp{site});
            end
            epsilon = abs(abs(nExp) - tc.LENGTH * (tc.HILBY - 1));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testNonHermState(tc)
            dmpo = tc.dmpo;
            midSite = ceil(tc.LENGTH/2);
            op = tc.idOp;

            op(:, :, midSite) = zeros(tc.HILBY);
            for state = 1 : 1 : (tc.HILBY - 1)
                op(state, state + 1, midSite) = sqrt(state);
            end
            exp0 = DMPOExp(dmpo, op);

            % tamper with the DMPO slightly to ensure different results
            % if the indices are the wrong way round
            dmpo{midSite}(:, :, 2, 1) = -1 * dmpo{midSite}(:, :, 2, 1);

            % the expectation value *should* now be different
            exp1 = DMPOExp(dmpo, op);
            tc.assertNotEqual(exp1, exp0);
        end
    end
end
