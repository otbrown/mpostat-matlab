% ConvTest.m
% function which checks for convergence of a set of values below some
% threshold -- uses the mean of the absolute difference between successive
% values
% Oliver Thomson Brown
% 2016-08-30
%
% [ convFlag, convergence ] = ConvTest(data, THRESHOLD)
%
% RETURN
% convFlag: 	bool, true if convergence is below threshold
% convergence: 	double, the mean absolute change between successive values
%				in data
%
% INPUT
% data: 		(complex) double array, the set of values to be tested for
%				convergence
% THRESHOLD: 	double, convergence threshold

function [ convFlag, convergence ] = ConvTest(data, THRESHOLD)
	convFlag = false;
	convergence = NaN;

	if ~isnan(data)
		convergence = mean(abs(diff(data)));
		if convergence < THRESHOLD
			convFlag = true;
		end
	end
end
