function [data, wires, x, f] = loadVcdData_moduled(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, overwrite, isbin2dec)
pb = CmdLineProgressBar('Reading VCD files');
filenameout = 'mode';
for m=1:length(Mvals)
filenameout=[filenameout num2str(Mvals(m))];  %#ok<*AGROW>
end
filenameout = [filenameout '_y_'];
for m=1:length(Yvals)
filenameout=[filenameout num2str(Yvals(m))]; %#ok<*NASGU>
end
filenameout = [filenameout '_x_' num2str(NparamsX) 'vals_z_' num2str(NparamsZ) 'vals_moduled.mat'];

if ~exist(filenameout,'file') || overwrite 
    
    inds2params = randperm(Nparams_all);
    Xvec = sort(uniqueX(inds2params(1:NparamsX)),'ascend');
    Zvec = sort(uniqueZ(inds2params(1:NparamsZ)),'ascend');  
    n=1;
    for Y=Yvals
        for M = Mvals
            for X=Xvec
                for Z=Zvec
                    filenameIn = ['X' num2str(X) '_Y' num2str(Y) '_Z' num2str(Z) '_MODE' num2str(M) '.vcd'];
                    if ~exist(fullfile(datapth, filenameIn), 'file')
                        continue;
                    end
                    pb.print(n, length(Yvals)*length(Mvals)*NparamsX*NparamsZ);
                    [data.samps{n}, wires{n}] = parse_vcd_moduled(fullfile(datapth, filenameIn),0, inf, isbin2dec);
                    wires{n}.name = zeros(length(wires{n}.module),1);
                   
                    data.X(n) = X;
                    data.Y(n) = Y;
                    data.Z(n) = Z;
                    data.M(n) = M;
                    n = n + 1;
                end
            end
        end
    end

    save(filenameout,'data','wires');
else
    load(filenameout);
end
x = [];f=zeros(length(data.samps), 1);
for n=1:length(data.samps)
    if length(data.samps{n})~=441
        f(n)=0;
    else
            x = cat(3, x, data.samps{n}(:,4:end));
        f(n)=1;
    end
end