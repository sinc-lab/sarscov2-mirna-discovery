function fold_sequences( nworks )
display('Folding...');
tic

parfor i=1:nworks
    system(['RNAfold --noPS < toerase_windows_' num2str(i) '.fasta > toerase_windows_' num2str(i) '.fold']);
end
t=toc;
display(['   - Elapsed time: ' num2str(t/60) ' min']);

end

