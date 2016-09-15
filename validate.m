% validate.m
% runs both unit testing and the validation tests
% Oliver Thomson Brown
% 2016-09-08

% set up result area
resultDir = 'testing/results';
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

utilResult = runtests('testing/utilTests');
coreResult = runtests('testing/coreTests');
topResult = runtests('testing/topTests');
valResult = runtests('testing/validationTests');

failCount = sum([utilResult.Failed]) + sum([coreResult.Failed]) ...
            + sum([topResult.Failed]) + sum([valResult.Failed]);

if failCount
    dateStamp = datestr(datetime('now'), 30);
    savePath = sprintf('%s/%s_VALIDATE.mat', resultDir, dateStamp);
    save(savePath);
end
