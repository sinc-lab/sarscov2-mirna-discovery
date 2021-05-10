function [outUnlabeled, outMiRNA, error_log] = HExtractor(inputFile, configFile, annotatedMiRNAs, otherRNA)
% Inizialization
clear global *;
rng('default');

wp = mfilename('fullpath');
wp = wp(1:end-11);
cd(wp);

if inputFile(1) ~= '/'
    inputFile = [wp filesep 'input' filesep inputFile];
end
if configFile(1) ~= '/'
    configFile = [wp filesep 'config' filesep  configFile];
end
if nargin>2 && ~isempty(annotatedMiRNAs) && annotatedMiRNAs(1) ~= '/'
    annotatedMiRNAs = [wp filesep 'input' filesep annotatedMiRNAs];
end
if nargin>3 && ~isempty(otherRNA) && otherRNA(1) ~= '/'
    otherRNA = [wp filesep 'input' filesep  otherRNA];
end

addpath(genpath(wp));
dp = [wp '/lib/YAMLMatlab/external/snakeyaml-1.13.jar'];
if not(ismember(dp, javaclasspath ('-dynamic')))
    javaaddpath(dp); % javaaddpath clears global variables...!?
end
addpath(genpath('src'));

% Opening files
config = ReadYaml(configFile);

filter_files = {};
if nargin > 2 && ~isempty(annotatedMiRNAs)
    filter_files = {annotatedMiRNAs};
end
if nargin > 3 && ~isempty(otherRNA)
    filter_files = [filter_files, {otherRNA}];
end

% Creating needed threads
if ~isfield(config, 'nthreads')
    config.nthreads = 1;
end
if ~isfield(config, 'nworks')
    config.nworks = 1;
end

config.nworks = min(config.nthreads, config.nworks);
if config.nthreads > 1
    ppool_open( config.nthreads )
end


numItems = split_fasta(inputFile);

go_tmp( wp );

status = zeros(1,numItems);
error = cell(1,numItems);
for iseq=1:numItems
    sequence = fastaread([inputFile, '.', num2str(iseq)]);
    sequence.Header = config.name_seq;
    if ~exist(['part_' sequence.Header '.fasta'], 'file')
        if sum(sequence.Sequence ~= 'N') < config.min_valid_nucleotides
            error{iseq} = ['Skipping sequence ' sequence.Header ' due to lack of valid nucleotides'];
            continue;
        end
        try
            delete('toerase_*');
            window_sequence(sequence.Sequence, sequence.Header, config.nworks, ...
               config.window_step, config.window_size, config.min_seq_length);
            
            fold_sequences( config.nworks );
            extract_stems( config.nworks, config.min_seq_length, ...
               config.min_bp, config.margin_bp, config.only_sloop);
            unique_sequences(config.nworks, ['part_' sequence.Header '.fasta']);
        catch err
            error{iseq} = ['Sequence ' sequence.Header ': ' err.message '\n'];
        end
    end
    delete([inputFile, '.', num2str(iseq)])
end
filter_result = exec_blast(filter_files, config.blast_evalue);
final_files = separate_sequences(filter_result, config.identity_threshold, [wp filesep 'out'], config.name_out);

outUnlabeled = final_files{1};
if length(final_files) == 2
    outMiRNA = final_files{2};
else
    outMiRNA = '';
end

error_log = [error{:}];
return_tmp();
end
