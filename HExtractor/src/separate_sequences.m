function final_files = separate_sequences( filter_result, identity_threshold, final_path, name_out)
display('Separating sequences according to the filter files...');
tic

fid = fopen('all_stems.fasta','rt');
blockSize = 2^24;

hNdx = [];
pos = 0;
count = blockSize;
while count == blockSize
    [str, count] = fread(fid,blockSize,'*char');
    str = str';
    hNdx = [hNdx (pos + strfind(str,'>'))];
    pos = ftell(fid);
end
frewind(fid);

numItems = numel(hNdx);
headers = cell(size(hNdx));

for count = 1:numItems
    fseek(fid,hNdx(count),-1);
    headers{count} = fgetl(fid);
end
fclose(fid);

label = ones(numItems,1);
for i=1:length(filter_result)
    d=table2cell(readtable(filter_result{i}, 'Delimiter',',','ReadVariableNames',false));
    if ~isempty(d)
        query=d(:,2);
        identity=[d{:,3}];
        query=query(identity>identity_threshold);
        query=sort(query);
        label(ismember(headers, query))=i+1;
    end
end

taux = clock();
timestamp = [num2str(taux(1)) '-' num2str(taux(2)) '-' num2str(taux(3))...
    '-' num2str(taux(4)) num2str(taux(5)) num2str(ceil(taux(6)))];

final_files = cell(1,length(filter_result));
final_ids = cell(1,length(filter_result));
%final_files{1} = [final_path filesep timestamp 'filtered_unlabeled.fasta'];
final_files{1} = [final_path filesep name_out '.fasta'];
final_ids{1} =  fopen(final_files{1},'w');

display(['   -  unlabeled: ' num2str(sum(label==1))]);
for i=2:(length(filter_result)+1)
    filename = filter_result{i-1}(1:(end-4));
    if i < 3
        final_files{i} = [final_path filesep name_out '_filtered-' filename '.fasta'];        
        final_ids{i} =  fopen(final_files{i},'w');
    end
    display(['   -  ' filename ': ' num2str(sum(label==i))]);
end

fid = fopen('all_stems.fasta','r');

for i = 1:numItems
    chunkStart = hNdx(i)-1;
    if i < numItems
        chunkEnd = hNdx(i+1)-1;
        chunkSize = chunkEnd - chunkStart;
    else
        chunkSize = Inf;
    end
    
    fseek(fid,chunkStart,-1);
    chunk = fread(fid,chunkSize,'*char');
    if label(i) <= length(final_ids)
        fwrite(final_ids{label(i)}, chunk, 'char');
    end
end
fclose(fid);

for i=1:length(final_ids)
    fclose(final_ids{i});
end

display(['   - Elapsed time: ' num2str(toc) ' sec']);
end
