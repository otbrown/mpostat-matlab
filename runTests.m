% runTests.m
% platform independent script to run test-suite on mpostat-matlab and save
% results, eases remote deployment by creating a single call for the suite
% and separating results
% Oliver Thomson Brown
% 2016-08-30

% set up result area
resultDir = 'testing/results';
mkdir(resultDir);
dateStamp = datestr(datetime('now'), 30);
savePath = sprintf('%s/%s_TEST.mat', resultDir, dateStamp);

result = runtests('testing', 'IncludeSubFolders', true);

save(savePath);
