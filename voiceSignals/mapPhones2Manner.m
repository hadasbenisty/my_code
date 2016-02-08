function [manner_str, manner_num] = mapPhones2Manner(phones)


manner_str = cell(size(phones));
manner_num=zeros(length(phones), 1);
% vowels 1
for n=1:length(phones)
    switch phones{n}
        % vowels
        case {'aa', 'ae',  'ao',  'ax',  'axh',...
         'axr',  'ay',  'iy',  'ih',  'eh',...
         'ey',   'oy',  'ow',  'uh',...
         'uw',  'ux',  'er',  'ix'}
        manner_str{n} = 'v';
        manner_num(n) = 1;
        case {'b',  'd',  'g',  'p',  't',  'k',  'dx',  'q'}
            manner_str{n} = 's';
        % Affricates+Fricatives
        case {'jh',  'f', 'ch', 's',  'sh',  'z',  'zh',...
         'th',  'v',  'dh'}
        manner_str{n} = 'f';
        manner_num(n) = 2;
        
        % Nasals 3
        case  {'m',  'n',  'ng',  'em',  'en',  'eng',  'nx'}
            manner_str{n}= 'n';
            manner_num(n) = 3;
        %   Semivowels and Glides 4
        case {'l',  'r',  'w',  'hh',  'hv',  'el','y'}
            manner_str{n} = 'g';
            manner_num(n) = 4;
        % sil and close 5
        case {'bcl', 'dcl',  'gcl',  'pcl',  'tcl',  'kcl',...
         'pau',  'epi',  'h#',  'ssil'}
        manner_str{n} = 'c';
        manner_num(n) = 5;
    end
end
