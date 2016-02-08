function [wordsSamples, wordsStr] = readWrdFile(filename, x)


fid = fopen(filename);
C = textscan(fid, '%d %d %s');
fclose(fid);

begin_time = C{1};
end_time = C{2};
wordsStr = C{3};






wordsSamples = cell(length(C{1}), 1);
for T = 1:length(C{1})
   wordsSamples{T} = x(max(begin_time(T), 1):min(end_time(T), length(x)));
end

