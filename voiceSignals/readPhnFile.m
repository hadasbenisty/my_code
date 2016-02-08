function phones = readPhnFile(filename, nfft, hop)


fid = fopen(filename);
C = textscan(fid, '%d %d %s');
fclose(fid);

begin_time = round(C{1}/nfft/hop)+1;
end_time = round(C{2}/nfft/hop);

phones = cell(end_time(end),1);
for n=1:length(end_time)
    for m=begin_time(n):end_time(n)
        phones{m} = C{3}{n};
    end
end