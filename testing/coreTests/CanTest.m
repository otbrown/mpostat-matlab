% CanTest.m
% Oliver Thomson Brown
% 2016-03-25

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev', 'IncludingSubfolders', true)}) CanTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        HILBY = 2;
        LENGTH = 5;
        COMPRESS = 0;
        dmpo;
    end

    methods (TestClassSetup)
        function MethodSetup(tc)
            tc.dmpo = ZDMPO(tc.HILBY, tc.LENGTH, tc.COMPRESS);
        end
    end

    methods (Test)
        function testTracePreservation(tc)
            ldmpo = Can(tc.dmpo, [1 : 1 : (tc.LENGTH-1)], 'L');
            rdmpo = Can(tc.dmpo, [tc.LENGTH : -1 : 2], 'R');
            tr = DMPOTrace(tc.dmpo);
            trL = DMPOTrace(ldmpo);
            trR = DMPOTrace(rdmpo);
            tc.assertLessThan(abs(tr - trL), tc.absTol);
            tc.assertLessThan(abs(tr - trR), tc.absTol);
        end

        function testCanBadDirection(tc)
            route = 1 : 1 : (tc.LENGTH-1);
            badDirection = 'm';
            tc.assertError(@()Can(tc.dmpo, route, badDirection), ...
                                    'Can:BadDirection');
        end

        function testLCanBadRoute(tc)
            direction = 'L';
            badSite = tc.LENGTH;
            badEndRoute = 1 : 1 : tc.LENGTH;
            wrongWayRoute = (tc.LENGTH-1) : -1 : 2;
            tc.assertError(@()Can(tc.dmpo, badSite, direction), ...
                                    'Can:LCan:BadRoute');
            tc.assertError(@()Can(tc.dmpo, badEndRoute, direction), ...
                                    'Can:LCan:BadRoute');
            tc.assertError(@()Can(tc.dmpo, wrongWayRoute, direction), ...
                                    'Can:LCan:BadRoute');
        end

        function testRCanBadRoute(tc)
            direction = 'R';
            badSite = 1;
            badEndRoute = tc.LENGTH : -1 : 1;
            wrongWayRoute = 2 : 1 : (tc.LENGTH - 1);
            tc.assertError(@()Can(tc.dmpo, badSite, direction), ...
                                    'Can:RCan:BadRoute');
            tc.assertError(@()Can(tc.dmpo, badEndRoute, direction), ...
                                    'Can:RCan:BadRoute');
            tc.assertError(@()Can(tc.dmpo, wrongWayRoute, direction), ...
                                    'Can:RCan:BadRoute');
        end
    end
end
