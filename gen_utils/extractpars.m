function pars=extractpars(vars,default)
% function pars=extractpars(vars,default);
%
%

if(nargin<2)
    default=[];
end

pars=default;
if(length(vars)==1)
    p=vars{1};
    s=fieldnames(p);
    for i=1:length(s)
        eval(['pars.' s{i} '=p.' s{i} ';']);
    end 
else
    for i=1:2:length(vars)
        if(ischar(vars{i}))
            if(i+1>length(vars)) error(sprintf('Parameter %s has no value\n',vars{i}));else val=vars{i+1};end;
            if(ischar(val))
                pars.(vars{i}) = val; %eval(['pars.' vars{i} '=''' val ''';']);
            else
                for j=1:length(val)
                    eval(['pars.' vars{i} '(' num2str(j) ')=' sprintf('%i',val(j)) ';']);
                end
            end
        end
    end 
end

