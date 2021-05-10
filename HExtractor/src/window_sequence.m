function window_sequence( seq, seqName, nworks, step, win, min_length )
display('- Windowing...');
tic

seq = strrep(lower(seq), 't', 'u');
sic = seq;
sic(seq=='a')='u';
sic(seq=='u')='a';
sic(seq=='c')='g';
sic(seq=='g')='c';

s_add = round(2*length(seq)/step);
to_add= -ones(s_add,2);
n_add = 0;
for i=1:step:(length(seq)-win)
    [ibeg, iend] = regexp(seq(i:i+win), '[acgu]+');
    for j=1:length(ibeg)
        if iend(j)-ibeg(j) > min_length
            n_add = n_add+1;
            if n_add<=s_add
                to_add(n_add,:)= [ibeg(j) iend(j)]-1+i;
            else
                to_add = [to_add; [ibeg(j) iend(j)]-1+i];
            end
        end
    end
end
to_add=to_add(1:n_add,:);
to_add=unique(to_add,'rows');
n_add = 2*size(to_add,1);

for i=1:nworks
    fid=fopen(['toerase_windows_' num2str(i) '.fasta'], 'w');
    
    fseq = 1+round((i-1)*size(to_add,1)/nworks);
    lseq = round(i*size(to_add,1)/nworks);
    for j=fseq:lseq
        fprintf(fid, '>%s_%i-%i\n%s\n', seqName, to_add(j,1), to_add(j,2), ...
            seq(to_add(j,1):to_add(j,2)));
    end
    
    for j=fseq:lseq
        fprintf(fid, '>%s_%i-%i_InvCom\n%s\n', seqName, to_add(j,2), to_add(j,1), ...
            sic(to_add(j,2):-1:to_add(j,1)));
    end
    fclose(fid);
end

t=toc;
display(['   -> Number of sequences: ' num2str(n_add)]);
display(['   -> Elapsed time: ' num2str(t) 'sec']);
end

