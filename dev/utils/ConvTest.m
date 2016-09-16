% ConvTest.m
% function which checks if the supplied values have converged below threshold
% Oliver Thomson Brown
% 2016-08-30
%
% [RETURN]
% convFlag	: test result, convFlag == 1.0 if successful and convFlag == 0.0 if
%             test is failed
%
% [INPUTS]
% eigTracker	: complex array, contains the eigenvalue of the effective
%                 Liouvillian solved for each site update
% threshold	    : double, convergence threshold

function [ convFlag, convergence ] = ConvTest(eigTracker, threshold)
	convFlag = false;
	convergence = NaN;

	if ~isnan(eigTracker)
		convergence = mean(abs(diff(eigTracker)));
		if convergence < threshold
			convFlag = true;
		end
	end
end
