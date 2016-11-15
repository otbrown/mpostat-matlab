% profiling.m
% platform independent script to profile variational stationary state finder
% use the command:
%   profview(0, profinfo)
% to view the profile results created in this script
% Oliver Thomson Brown
% 2016-09-09

addpath(genpath('dev'));
addpath('external/primme/Matlab');

% set up result area
resultDir = 'profiling';
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

% set up system to be solved
HILBY = 3;
LENGTH = 4;
COMPRESS = 75;
THRESHOLD = 1E-5;
variant = 'direct';

dateStamp = datestr(datetime('now'), 30);
savePath = sprintf('%s/%s_%sPROFILE.mat', resultDir, dateStamp, variant);

dmpoInit = MixDMPO(HILBY, LENGTH, COMPRESS);

% build Liouvillian mpo
% parameters
detuning01 = 0;
detuning02 = 0;
hop = 0.5;
diss21 = 1;
diss10 = 0.1;

% local operators
I = eye(HILBY);
a = zeros(HILBY);
for i = 1 : 1 : (HILBY - 1)
    a(i, i + 1) = sqrt(i);
end
a10 = [0, 1, 0; 0, 0, 0; 0, 0, 0];
a21 = [0, 0, 0; 0, 0, sqrt(2); 0, 0, 0];

% mpo building blocks
ident = kron(I, I);
Hloc = [0, 0, 0; 0, detuning01, 0; 0, 0, detuning02];
dissipator = 0.5 * diss10 * (kron(2 * conj(a10), a10) ...
             - kron(I, ctranspose(a10) * a10) ...
             - kron(transpose(a10) * conj(a10), I)) ...
             + ...
             0.5 * diss21 * (kron(2 * conj(a21), a21) ...
             - kron(I, ctranspose(a21) * a21) ...
             - kron(transpose(a21) * conj(a21), I));

% put liouvillian mpo together
mpoLoc = kron(I, -1i*Hloc) + kron(1i*Hloc, I) + dissipator;

lmpo = zeros(HILBY, HILBY, HILBY, HILBY, 6, 6);
lmpo(:, :, :, :, 6, 1) = reshape(mpoLoc, ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 1, 1) = reshape(ident, ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 6, 6) = reshape(ident, ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 2, 1) = reshape(kron(I, a), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 3, 1) = reshape(kron(a, I), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 4, 1) = reshape(kron(I, ctranspose(a)), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 5, 1) = reshape(kron(ctranspose(a), I), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 6, 2) = reshape(kron(I, 1i*hop*ctranspose(a)), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 6, 3) = reshape(kron(-1i*hop*ctranspose(a), I), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 6, 4) = reshape(kron(I, 1i*hop*a), ...
                         [HILBY, HILBY, HILBY, HILBY]);
lmpo(:, :, :, :, 6, 5) = reshape(kron(-1i*hop*a, I), ...
                         [HILBY, HILBY, HILBY, HILBY]);

% build mpo cell
mpo = cell(LENGTH, 1);
mpo{1} = lmpo(:, :, :, :, 6, :);
mpo{LENGTH} = lmpo(:, :, :, :, :, 1);
for site = 2 : 1 : (LENGTH - 1)
    mpo{site} = lmpo;
end

profile on;

% solve using Stationary
[dmpoStat, eigTrack] = Stationary(dmpoInit, mpo, THRESHOLD, variant);

profinfo = profile('info');

% a couple of basic checks that the solution is correct
% this system should end up in the all zero state, so we'll check trace and
% the first element of the density matrix
tr = DMPOTrace(dmpoStat);
rho0 = 1;
for site = 1 : 1 : LENGTH
    rho0 = rho0 * dmpoStat{site}(:, :, 1, 1);
end

save(savePath);
