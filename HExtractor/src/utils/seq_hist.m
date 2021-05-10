function seq_hist( file, varargin )
minima=50;
maxima=200;
if(nargin == 3)
    minima=varargin{1};
    maxima=varargin{2};
end

if iscell(file)
    seq = fastaread(file{1});
    for i=2:length(file);
        seq = [seq; fastaread(file{2})];
    end
else
    seq = fastaread(file);
end

len = zeros(length(seq),1);
for i=1:length(seq)
    len(i) = length(seq(i).Sequence);
end
len = len(len<maxima);

hist(len, minima:5:maxima);
end
