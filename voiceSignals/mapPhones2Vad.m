function [vad_str, vad_num] = mapPhones2Vad(phones)


vad_str = cell(size(phones));
vad_num=zeros(length(phones), 1);
% vowels 1
for n=1:length(phones)
    switch phones{n}
       
        case {'bcl', 'dcl',  'gcl',  'pcl',  'tcl',  'kcl',...
         'pau',  'epi',  'h#',  'ssil'}
        vad_str{n} = 'sil';
        vad_num(n) = 0;
        otherwise
             vad_str{n} = 'sp';
        vad_num(n) = 1;
    end
end
