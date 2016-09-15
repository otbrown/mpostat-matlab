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
val4SResult = runtests('testing/validationTests', 'Name', 'ValidationTest/FourSiteExact');
valMFS15Result = runtests('testing/validationTests', 'Name', 'ValidationTest/MFS15')

failCount = sum([utilResult.Failed]) + sum([coreResult.Failed]) ...
            + sum([topResult.Failed]) + sum([val4SResult.Failed]) + sum([valMFS15Result.Failed]);
if failCount
    dateStamp = datestr(datetime('now'), 30);
    savePath = sprintf('%s/%s_VALIDATE.mat', resultDir, dateStamp);
    save(savePath);
end
