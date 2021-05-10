function extract_stems( nworks, min_length, min_bp, margin, only_sloop )
display('Extracting stems...');
tic

parfor it=1:nworks
    fid=fopen(['toerase_stems_' num2str(it) '.fasta'], 'w');
    in_seqs = fastaread(['toerase_windows_' num2str(it) '.fold']);
    for i=1:length(in_seqs)
        [seq,str,par] = parseRNAfoldOut(in_seqs(i).Sequence);
        
        if sum(str~='.') < 5
            continue
        end
        
        if(only_sloop)
            cut = locateCutPositions(str,par);
            [ibeg, iend] = getIndexs(cut, str, only_sloop);
            for j=1:length(ibeg)
                if (iend(j)-ibeg(j)) > min_length && sum(str(ibeg(j):iend(j))=='(') >= min_bp
                    sseq=seq((ibeg(j):iend(j)));
                    %sstr=str(ibeg(j):iend(j));
                    %spar=par(ibeg(j):iend(j)) - ibeg(j) + 1;
                    %sseq=trimSequence(sseq,sstr,spar,min_bp,margin);
                    sname=in_seqs(i).Header
                    sname(sname==' ')='_';
                    fprintf(fid, '>%s_stem-%i-%i\n%s\n', sname, ibeg(j), iend(j), sseq);
                end
            end
        else
            if length(seq) > min_length && sum(str)~='.' > 2*min_bp
                fprintf(fid, '>%s-%i\n%s\n', in_seqs(i).Header, i, seq);
            end
        end
    end
    fclose(fid);
end
t=toc;
display(['   - Elapsed time: ' num2str(t) ' sec']);
end

function [seq] = trimSequence(seq,str,par,minbp,margin)
if find(str=='(', 1, 'last') > find(str==')', 1, 'first')
    return;
end

mid = round((regexp(str, '\(\.*\)', 'start')+regexp(str, '\(\.*\)', 'end'))/2);
flen=zeros(1,length(seq));
dist=0;
for i=mid:length(seq)
    flen(i)=dist;
    if( str(i) ~= '.' )
        dist=dist+1;
    end
end
dist=0;
for i=mid:-1:1
    flen(i)=dist;
    if( str(i) ~= '.' )
        dist=dist+1;
    end
end
p1 = (flen-minbp-margin) > 0; % aca resta un vector de distancias por el margen y el min de bp de largo

db = regexp(str, '\.*', 'start');
de = regexp(str, '\.*', 'end');
flen=zeros(1,length(seq));
for i=1:length(de)
    flen(db(i):de(i)) = de(i)-db(i); % busca los \.* y guarda sus longitudes a lo largo de la longitud del mismo
end 
p2 = flen;

flen = -abs((1:length(p2))-mid);
p3 = (1-flen/min(flen)).^2; %porque divide por el minimo?? deberia hacerlo por el maximo

p = p1 .* p2 .* p3;
if(any(p>0))
    i = find(max(p)==p);
    sbeg = i(min(abs(i-mid))==abs(i-mid));
    sbeg = sbeg(1);
    send = max(par(sbeg+1), par(sbeg-1));
    if sbeg<send
        seq = seq(sbeg:send);
    else
        seq = seq(send:sbeg);
    end
end
end

function [seq, str, par] = parseRNAfoldOut(string)
seq=char(regexpi(string, '[acgut]+', 'match', 'once'));
str=char(regexpi(string, '[(.)]+', 'match', 'once'));
par=-ones(1,length(str));
stack=zeros(1,length(str));
ss=0;
for j=1:length(par)
    if str(j)=='('
        ss=ss+1;
        stack(ss)=j;
    elseif str(j)==')'
        par(stack(ss))=j;
        par(j)=stack(ss);
        ss=ss-1;
    end
end
end

function cuts = locateCutPositions(str, par)
marks = [find(str=='(',1,'first'), regexpi(str, '\)\.*\(', 'start'), ...
    find(str==')',1,'last'),  regexpi(str, '\)\.*\(', 'end')];
cuts=false(1,length(str));
cuts(marks)=true;
cuts(par(marks))=true;
end

function [ibeg, iend] = getIndexs(cut, str, only_sloop)
ibeg = -ones(1,sum(cut)/2);
iend = ibeg;
nadded = 0;
if only_sloop
    last_cut=-1;
    for i=find(cut)
        if str(i)=='('
            last_cut=i;
        elseif str(i)==')' && last_cut>0
            nadded=nadded+1;            
            ibeg(nadded)=last_cut-1;
            
            while ibeg(nadded)>0 && str(ibeg(nadded))=='.'
                ibeg(nadded)=ibeg(nadded)-1;
            end
            
            ibeg(nadded)=ibeg(nadded)+1;
            
            iend(nadded)=i+1;
            
            while iend(nadded)<=length(str) && str(iend(nadded))=='.'
                iend(nadded)=iend(nadded)+1;
            end
            
            iend(nadded)=iend(nadded)-1;
            if ((iend(nadded)+10)<=length(str)) % nc
                iend(nadded)=iend(nadded)+10;   % nc
            end                                 % nc
                
            last_cut=-1;
        end
    end
else
    stack=zeros(1,length(str));
    ss=0;
    for i=find(cut)
        if str(i)=='('
            ss=ss+1;
            stack(ss)=i;
        elseif str(i)==')'
            nadded=nadded+1;
            ibeg(nadded)=max(stack(ss)-5, 1);
            iend(nadded)=min(i+5,length(str));
            ss=ss-1;
        end
    end
end
end
