function unique_sequences( nworks, ofilename )
display('Filtering repeated sequences...');
tic

final_sequences = [];
nvector=zeros(1,nworks);
parfor i=1:nworks
    seqs = fastaread(['toerase_stems_' num2str(i) '.fasta']);
    nvector(i) = length(seqs);
    
    [~, idx] = unique({seqs.Sequence});
    idx=sort(idx);

    final_sequences = [final_sequences; seqs(idx)];
end
n1 = sum(nvector);

[~, idx] = unique({final_sequences.Sequence});
idx=sort(idx);
final_sequences=final_sequences(idx);
n2=length(final_sequences);

uni = true(1,length(final_sequences));
for i=11:length(uni)-11
    if ~uni(i)
        continue;
    end
    atom = final_sequences(i).Sequence;
    atom = atom(floor(length(atom)/2-10):ceil(length(atom)/2+10));
    for j=i-10:i+10
        if i==j || ~uni(j)
            continue
        end
        rgx = regexp(final_sequences(j).Sequence, atom, 'once');
        if ~isempty(rgx) && rgx > 10 && rgx < length(final_sequences(j).Sequence)-10
            if length(final_sequences(i).Sequence) > length(final_sequences(j).Sequence)
                uni(j)=false;
            else
                uni(i)=false;
            end
        end
    end
    
end
final_sequences=final_sequences(uni);
n3=length(final_sequences);

if exist(ofilename, 'file')
    delete(ofilename)
end
fastawrite(ofilename, final_sequences);

display(['   - Deleted sequences: ' num2str(n1) ' -> ' num2str(n2) ' -> ' num2str(n3)]);
display(['   - Elapsed time: ' num2str(toc) ' sec']);
