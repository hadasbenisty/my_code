function [data, wires] = parse_vcd_moduled(filename, verbose, numvars, isbin2dec)

fid = fopen(filename,'r');
if nargin < 3
    numvars = inf;
    if nargin < 2
        verbose = 0;
    end
end
MAX_LINES = 1e6;
MAX_VARS = 1e3;
MAX_TIME = 1e5;
data=zeros(MAX_VARS, MAX_TIME);
state = 0;% 0 - var names, 1 - wires val, 2 - error, 3 - finish
wires.name = [];
wires.width = [];
wires.label = [];
wires.module=[];
wires.mdl_name = [];
wires.widthind = [];
curr_t = 1;
currmdle=[];
wires.time = [];
for i_line = 1:MAX_LINES
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if verbose > 0
        disp(tline);
    end
    switch state
        case 0
            if ~isempty(strfind(tline,'$enddefinitions $end'))
                state = 1;
                tline = fgetl(fid);
                tline = fgetl(fid);
                tline = fgetl(fid);
                if ~isinf(numvars) && length(wires.name) ~= numvars
                    error('numvars is not equal to the amount of read wires');
                end
                data = data(1:length(wires.name),:);
                continue;
            end
            
        if ~isempty(strfind(tline, '$upscope $end'))            
            currmdle = currmdle(1:end-1);
            continue;
        end
        if ~isempty(strfind(tline, '$scope module'))
            c = textscan(tline,'$scope module %s $end');
            currmdle{end+1} = c{1}{1};
            continue;
        end
        if isempty(strfind(tline,'$var'))
            continue;
        end
        c=textscan(tline,'%s');
        if isbin2dec
        wires.width(end+1) = str2double(c{1}{3});
        wires.label{end+1} = c{1}{4};
        wires.name{end+1} = c{1}{5};
        wires.module{end+1} = currmdle;
        wires.mdl_name{end+1} = [wires.module{end} ', ' wires.name{end}];
        else
            for w = 1:str2double(c{1}{3})
                wires.width(end+1) = 1;
                wires.widthind(end+1) = w;
                wires.label{end+1} = [c{1}{4} ];
                wires.name{end+1} = c{1}{5} ;
                wires.module{end+1} = currmdle;
                wires.mdl_name{end+1} = [wires.module{end} ', ' wires.name{end}];
            end
        end
        case 1
            if strcmp(tline(1), '#') && isnumeric(str2double(tline(2:end)))
                wires.time(end+1) = str2double(tline(2:end));
                curr_t=curr_t+1;
                continue;
            end
            c=textscan(tline,'%s');
            if strcmp(tline,'$dumpon');
                continue;
            end
            if strcmp(tline,'$end')
                curr_t=curr_t+1;
                continue;
            end
            if length(c{1}) ==1
                c1=c{1}{1};
                id = c1(2:end);
                valstr = c1(1);
            else
                id = c{1}{2};
                valstr = c{1}{1};
            end
            ind2id = find(strcmp(id,wires.label));
            if isempty(ind2id)
                error('Unidentified id');
            end
            if any(wires.widthind(ind2id)>1)
                for w=2:length(valstr)
                    data(ind2id(w-1), curr_t:end) = setvalFromVcdStr(valstr(w), isbin2dec);
                end
            else
                data(ind2id, curr_t:end) = setvalFromVcdStr(valstr, isbin2dec);
            end
            case 2
                case 3
end
end
fclose(fid);
data = data(:, 1:curr_t);
end


function numval = setvalFromVcdStr(vcdStr, isbin2dec)

switch vcdStr(1)
    case 'x'
        numval = nan;
    case {'b','B'} % binary
        inds = strfind(vcdStr(2:end),'z');
        vcdStr(inds+1) = '0';
        if isbin2dec
            if isempty(strfind(vcdStr(2:end),'x'))
                numval = bin2dec(vcdStr(2:end));
            else
                numval = nan;
            end
        else
            for w=2:length(vcdStr)
                numval(w-1) = my_str2num(vcdStr(w));
            end
        end
    otherwise
        if isbin2dec
            numval = str2double(vcdStr);
        else
            numval = my_str2num(vcdStr);
        end
end
end
function numval = my_str2num(vcSdtr)
switch vcSdtr
    case {'0', 'z'}
        numval= 0;
    case '1'
        numval = 1;
    case 'x'
        numval = inf;
    otherwise
        error('unrecongnized value');
end
end