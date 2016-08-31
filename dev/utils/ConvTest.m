
% ConvTest.m
% function which checks if the last sampleSize eigenvalues have been below
% threshold
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
% updCount      : double, the number of site updates performed so far
% sampleSize	: double, the number of results to check for convergence
% threshold	    : double, convergence threshold

function [ convFlag ] = ConvTest(eigTracker, updCount, sampleSize, threshold)
	convFlag = 0;

	if updCount <= sampleSize
        % total number of updates is still less than sampleSize
        % do nothing
	else
		sample = eigTracker(updCount - sampleSize : updCount) - ...
                 eigTracker(updCount);
		worstCase = max( abs(sample) );
		if worstCase < threshold
			convFlag = 1;
		end
	end
end
