function [voiced_str, voiced_num] = mapPhones2voicedUnvoices(phones)


voiced_str = cell(size(phones));
voiced_num=zeros(length(phones), 1);
% vowels 1
for n=1:length(phones)
    switch phones{n}
        % vowels
        case {'aa', 'ae',  'ao',  'ax',  'axh',...
         'axr',  'ay',  'iy',  'ih',  'eh',...
         'ey',   'oy',  'ow',  'uh',...
         'uw',  'ux',  'er',  'ix','b',  'd',  'g',  'p', ...
         't',  'k',  'dx',  'q','m',  'n',  'ng',  'em',  'en',  'eng',  'nx',...
         'l',  'r',  'w',  'hh',  'hv',  'el','y'}
        voiced_str{n} = 'v';
        voiced_num(n) = 1;
       
        case {'jh',  'f', 'ch', 's',  'sh',  'z',  'zh',...
         'th',  'v',  'dh', 'bcl', 'dcl',  'gcl',  'pcl',  'tcl',  'kcl',...
         'pau',  'epi',  'h#',  'ssil'}
        voiced_str{n} = 'u';
        voiced_num(n) = 0;
        
        
    end
end
