function mov = matrix2mov(data)

nframe = size(data, 3);
mov(1:nframe)=struct('cdata',[],'colormap',[]);
figure;
for k=1:nframe
imagesc(data(:, :, k));
mov(k)=getframe(gcf);
end