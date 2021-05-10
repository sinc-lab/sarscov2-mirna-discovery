function psize = ppool_size( )
%PPOOL_SIZE Summary of this function goes here
%   Detailed explanation goes here
v=version;
if str2num(v(1))<8
    if exist('matlabpool', 'file')==2 && matlabpool('size') > 0
        psize= matlabpool('size');
    else
        psize= 1;
    end
else
    if exist('gcp', 'file')==2
        pool = gcp('nocreate');
        if isempty(pool)
            psize=1;
        else
            psize= pool.NumWorkers;
        end
    else
        psize=1;
    end
end

