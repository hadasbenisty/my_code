function circ_params = initial_circ_files_analysis(circfile, datapth)
if ~exist(circfile,'file')
    n = 1;
    pb = CmdLineProgressBar('Reading VCD files');
    
    for Y=[9 58 255]
        for M = [0 1]
            for X=0:255
                for Z=0:255
                    filename = ['X' num2str(X) '_Y' num2str(Y) '_Z' num2str(Z) '_MODE' num2str(M) '.vcd'];
                    if ~exist(fullfile(datapth, filename), 'file')
                        continue;
                    end
                    pb.print(n, 60e3);
                    circ_params.X(n) = X;
                    circ_params.Y(n) = Y;
                    circ_params.Z(n) = Z;
                    circ_params.M(n) = M;
                    n = n + 1;
                end
            end
        end
    end
    save(circfile,'circ_params');
else
    load(circfile)
end
