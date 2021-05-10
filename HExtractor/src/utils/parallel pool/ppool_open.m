function ppool_open( nthreads )
%PPOOL_OPEN Summary of this function goes here
%   Detailed explanation goes here
v=version;
if str2num(v(1))<8
    if exist('matlabpool', 'file')==2 && matlabpool('size') < 1
        matlabpool OPEN nthreads;
    else
        error('Parallel toolbox not found');
    end
elseif exist('gcp', 'file')==2
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(nthreads);
    end
else
    error('Parallel toolbox not found');
end