function ppool_close( )
%PPOOL_CLOSE Summary of this function goes here
%   Detailed explanation goes here
v=version;
if str2num(v(1))<8
    if exist('matlabpool', 'file')==2 && matlabpool('size') > 0
        matlabpool CLOSE;
    end
elseif exist('gcp', 'file')==2
    pool = gcp('nocreate');
    if ~isempty(pool)
        delete(pool);
    end
end