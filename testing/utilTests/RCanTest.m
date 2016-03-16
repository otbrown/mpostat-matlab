% RCanTest.m
% Oliver Thomson Brown
% 2016-03-16

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) RCanTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        HILBY = 2;
        LENGTH = 7;
        COMPRESS = 100;
        SAMPLE_MAX = 100;
        dmpo;
        canDMPO;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 3, 4, 5};
        testLENGTH = {7, 5, 4, 5};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(tc, testHILBY, testLENGTH)
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
            tc.dmpo = DMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
            tc.canDMPO = RCan(tc.dmpo, [tc.LENGTH : -1 : 2]);
        end
    end

    methods (Test)
        % you might be thinking, "this looks a LOT like LCanTest". Yup.
        function testClass(tc)
            tc.fatalAssertClass(tc.canDMPO, 'cell');
        end

        function testSystemSize(tc)
            tc.fatalAssertSize(tc.canDMPO, size(tc.dmpo));
        end

        function testThrowBadRoute(tc)
            % this should be expanded at some point to check the route goes
            % through the system in the right direction
            badRoute = [tc.LENGTH : -1 : 1];
            tc.fatalAssertError(@()RCan(tc.dmpo, badRoute), 'RCan:BadRoute');
        end

        function testTensorShape(tc)
            % RCan may reduce virtual dimensions, but should still be consistent
            rowSz = 1;
            for site = 1 : 1 : (tc.LENGTH - 1)
                colSz = size(tc.canDMPO{site + 1}, 1);
                tc.fatalAssertSize(tc.canDMPO{site}, [rowSz, colSz, tc.HILBY, tc.HILBY]);
                rowSz = colSz;
            end
        end

        function testTracePreservation(tc)
            tr1 = DMPOTrace(tc.dmpo);
            tr2 = DMPOTrace(tc.canDMPO);
            epsilon = abs(abs(tr1) - abs(tr2));
            tc.assertLessThan(epsilon, tc.absTol);
        end

        function testCan(tc)
            % FIX ME
            %{
            for site = 1 : 1 : tc.LENGTH - 1
                [rowSz, colSz, ~, ~] = size(tc.canDMPO{site});
                M = reshape(tc.canDMPO{site}, [rowSz, colSz, tc.HILBY^2]);
                M = permute(M, [1, 3, 2]);
                M = reshape(M, [rowSz * tc.HILBY^2, colSz]);

                I = eye(colSz, colSz);
                epsilon = (ctranspose(M) * M) - I;
                tc.assertLessThan(epsilon, tc.absTol);
            end
            %}
        end

        function testElements(tc)
            % RCan should not alter the density matrix being represented
            SPACE = tc.HILBY^tc.LENGTH;
            sampleSz = min(floor(0.1 * SPACE^2), tc.SAMPLE_MAX);
            for testNum = 1 : 1 : sampleSz
                braState = randi([0, SPACE - 1]);
                ketState = randi([0, SPACE - 1]);
                braBits = FWBase(braState, tc.HILBY, tc.LENGTH) + 1;
                ketBits = FWBase(ketState, tc.HILBY, tc.LENGTH) + 1;

                coefft = 1;
                canCoefft = 1;
                for site = 1 : 1 : tc.LENGTH
                    bra = braBits(site);
                    ket = ketBits(site);
                    coefft = coefft * tc.dmpo{site}(:, :, bra, ket);
                    canCoefft = canCoefft * tc.canDMPO{site}(:, :, bra, ket);
                end

                epsilon = abs(abs(coefft) - abs(canCoefft));
                tc.assertLessThan(epsilon, tc.absTol);
            end
        end
    end
end
