function [new_phones_str, new_phones_num] = standardTIMITphones61to39(phones)

phoneLut = getFestvoxPhoneLut;

new_phones_str=phones;
new_phones_num = zeros(length(phones), 1);
for k=1:length(phones)
    if strcmp(phones{k}, 'ao')
        new_phones_str{k} = 'aa';
    elseif strcmp(phones{k}, 'ix')
        new_phones_str{k} = 'ih';
         elseif strcmp(phones{k}, 'dx')
        new_phones_str{k} = 't';
    elseif strcmp(phones{k}, 'el')
        new_phones_str{k} = 'l';
    elseif strcmp(phones{k}, 'em')
        new_phones_str{k} = 'm';
    elseif strcmp(phones{k}, 'zh')
        new_phones_str{k} = 'sh';
    elseif strcmp(phones{k}, 'ux')
        new_phones_str{k} = 'uw';
    elseif strcmp(phones{k}, 'ax') || strcmp(phones{k}, 'ax-h')
        new_phones_str{k} = 'ah';
    elseif strcmp(phones{k}, 'nx') || strcmp(phones{k}, 'en')
        new_phones_str{k} = 'n';
    elseif strcmp(phones{k}, 'bcl') || strcmp(phones{k}, 'pcl')||...
            strcmp(phones{k}, 'dcl')|| strcmp(phones{k}, 'tcl')||...
            strcmp(phones{k}, 'gcl')|| strcmp(phones{k}, 'kcl')|| ...
            strcmp(phones{k}, 'epi')|| strcmp(phones{k}, 'pau')|| ...
            strcmp(phones{k}, 'h#')
        new_phones_str{k} = 'pau';
    elseif strcmp(phones{k}, 'axr')
        new_phones_str{k} = 'er';
    elseif strcmp(phones{k}, 'hv')
        new_phones_str{k} = 'hh';
    elseif strcmp(phones{k}, 'eng')
        new_phones_str{k} = 'ng';
    elseif strcmp(phones{k}, 'q')
        new_phones_str{k} = 'k';
    end
    new_phones_num(k) = find(strcmp(phoneLut, new_phones_str{k}));
end

 
