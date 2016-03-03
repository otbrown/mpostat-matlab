% DMPOTraceTest.m
% Oliver Thomson Brown
% 2016-03-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../../dev')}) DMPOTraceTest < matlab.unittest.TestCase

    properties
        COMPRESS = 100;
        HILBY;
        LENGTH;
        dmpo;
        tr;
    end

    properties (MethodSetupParameter)
        testHILBY = {2, 2, 2, 3};
        testLENGTH = {5, 6, 7, 4};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function AssignParameters(tc, testHILBY, testLENGTH)
            % assign HILBY and LENGTH parameters to test case
            tc.HILBY = testHILBY;
            tc.LENGTH = testLENGTH;
        end

        function AllOnesDMPO(tc)
        	% COMPRESS == 0 means no compression
        	if tc.COMPRESS == 0
        		tc.COMPRESS = Inf;
        	end

            % allocate return
        	tc.dmpo = cell(tc.LENGTH, 1);

        	% first and last site
            tc.dmpo{1} = zeros(1, tc.HILBY^2, tc.HILBY, tc.HILBY);
            tc.dmpo{tc.LENGTH} = zeros(tc.HILBY^2, 1, tc.HILBY, tc.HILBY);
            for bra = 1 : 1 : tc.HILBY
                for ket = 1 : 1 : tc.HILBY
        	        tc.dmpo{1}(:, :, bra, ket) = eye(1, tc.HILBY^2);
        	        tc.dmpo{tc.LENGTH}(:, :, bra, ket) = eye(tc.HILBY^2, 1);
                end
            end

        	% and the rest
        	colSz = tc.HILBY^2;
        	for site = 2 : 1 : tc.LENGTH - 1
        		if site < ceil(tc.LENGTH / 2)
        			len = site;
        		else
        			len = tc.LENGTH - site;
        		end

        		rowSz = colSz;
        		colSz = min(tc.HILBY^(2*len), tc.COMPRESS);

        		tc.dmpo{site} = zeros(rowSz, colSz, tc.HILBY, tc.HILBY);
                for bra = 1 : 1 : tc.HILBY
                    for ket = 1 : 1 : tc.HILBY
                        tc.dmpo{site}(:, :, bra, ket) = eye(rowSz, colSz);
                    end
                end
        	end
        end

        function TraceCalc(tc)
            tc.tr = DMPOTrace(tc.dmpo);
        end
    end

    methods (Test)
        function testClass(tc)
            tc.fatalAssertClass(tc.tr, 'double');
        end

        function testReal(tc)
            % should be true in general, but will not be for highly entangled
            % highly compressed dmpos -- should work for the all ones density
            % matrices being used here however
            tc.fatalAssertEqual(imag(tc.tr), 0);
        end

        function testTrace(tc)
            % density matrix being built is just a 1 in every element -- as such
            % trace should equal HILBY^LENGTH
            SPACE = tc.HILBY^tc.LENGTH;
            tc.assertEqual(tc.tr, SPACE);
        end
    end
end
