function numItems = split_fasta(filename)
fid = fopen(filename,'rt');
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

fid = fopen(filename,'r');
for fileNum = 1:numItems
    outName = sprintf('%s.%d',filename,fileNum);
    fidOut = fopen(outName,'w');
    chunkStart = hNdx(((fileNum-1))+1)-1;
    if fileNum < numItems
        chunkEnd = hNdx(((fileNum))+1)-1;
        chunkSize = chunkEnd - chunkStart;
    else
        chunkSize = Inf;
    end
    fseek(fid,chunkStart,-1);
    chunk = fread(fid,chunkSize,'*char');
    fwrite(fidOut,chunk,'char');
    fclose(fidOut);
end
fclose(fid);
