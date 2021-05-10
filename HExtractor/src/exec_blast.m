function filter_result = exec_blast( filter_files, evalue )
display('Joining all stems and executing blast...');
tic

if ~exist('all_stems.fasta', 'file')
    system('cat part_* > all_stems.fasta');
    delete('part_*');
end
system('makeblastdb -max_file_sz 2GB -in all_stems.fasta -dbtype nucl > /dev/null');

filter_result = {};
for i=1:length(filter_files)
    filename = filter_files{i};
    uncompressedFile = {};
    if strcmp(filename(end-2:end), 'zip')
        uncompressedFile = unzip(filename);
        filename = uncompressedFile{1};
    end
    
    ini_filename = find(filter_files{i} == filesep, 1, 'last')+1;
    end_filename = find(filter_files{i} == '.', 1, 'last')-1;
    filter_result = [filter_result ...
        {[filter_files{i}(ini_filename:end_filename) '.csv']}];
    
    system(['blastn -query ' filename ' -db all_stems.fasta -evalue ' ...
        num2str(evalue) ' -dust no -strand plus -outfmt 10 -out ' filter_result{i}]);
    
    for j=1:length(uncompressedFile)
        delete(uncompressedFile{j})
    end
end

t=toc;
display(['   - Elapsed time: ' num2str(t/60) ' min']);
end
