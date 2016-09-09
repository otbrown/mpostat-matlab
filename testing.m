% testing.m
% platform independent script to run test-suite on mpostat-matlab and save
% results, eases remote deployment by creating a single call for the suite
% and separating results
% Oliver Thomson Brown
% 2016-08-30

% set up result area
resultDir = 'testing/results';
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

utilResult = runtests('testing/utilTests');
coreResult = runtests('testing/coreTests');
topResult = runtests('testing/topTests');

failCount = sum([utilResult.Failed]) + sum([coreResult.Failed]) ...
            + sum([topResult.Failed]);
if failCount
    dateStamp = datestr(datetime('now'), 30);
    savePath = sprintf('%s/%s_TESTING.mat', resultDir, dateStamp);
    save(savePath);
end
