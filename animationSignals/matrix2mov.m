function mov = matrix2mov(data, filename, ttl)
addpath(genpath('../../Matlab_Utils/'));
nframe = size(data, 3);
mov(1:nframe)=struct('cdata',[],'colormap',[]);
if nargin == 2
ttl = '';
end
figure;
for k=1:nframe
    imagesc(data(:, :, k));
title(ttl);
colormap gray;
    mov(k)=getframe(gcf);
end

if nargin > 1 && ~isempty(filename)
savemove2avi(filename, mov);
    
end