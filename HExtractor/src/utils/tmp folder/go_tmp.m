function go_tmp(script_path)
%GO_TMP move to a temporary directory
global miRNAfe_TMP_DIRECTORY;
global miRNAfe_OLD_DIRECTORY;

miRNAfe_OLD_DIRECTORY = pwd;
miRNAfe_TMP_DIRECTORY = [tempdir 'miRNAse-' num2str(int32(rand*2^31))];
while exist(miRNAfe_TMP_DIRECTORY, 'dir')
    miRNAfe_TMP_DIRECTORY = [tempdir 'miRNAse-' num2str(int32(rand*2^31))];
end

mkdir(miRNAfe_TMP_DIRECTORY);
cd(miRNAfe_TMP_DIRECTORY);
if nargin==1
    addpath(genpath(script_path));
end