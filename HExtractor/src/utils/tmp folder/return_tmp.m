function return_tmp(  )
%RETURN_TMP Summary of this function goes here
%   Detailed explanation goes here
global miRNAfe_TMP_DIRECTORY;
global miRNAfe_OLD_DIRECTORY;
if( isempty(miRNAfe_OLD_DIRECTORY) || isempty(miRNAfe_TMP_DIRECTORY))
    error('Can return to the pwd directory while already is in it');
end
%cd(miRNAfe_OLD_DIRECTORY);
% rmdir(miRNAfe_TMP_DIRECTORY,'s');

% clear global miRNAfe_TMP_DIRECTORY;
% clear global miRNAfe_OLD_DIRECTORY;
