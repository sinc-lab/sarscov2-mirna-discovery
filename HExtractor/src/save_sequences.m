function save_sequences( nworks, ofilename )
display('Saving sequences...');
tic

final_sequences = [];
for i=1:nworks
    final_sequences = [final_sequences; fastaread(['toerase_unique_' num2str(i) '.fasta'])];
end
[~, idx] = unique({final_sequences.Sequence});
if exist(ofilename, 'file')
    delete(ofilename)
end
fastawrite(ofilename, final_sequences(idx));

display(['   - Elapsed time: ' num2str(toc) ' sec']);
