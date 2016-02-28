function mov = matrix2mov(data, filename)

nframe = size(data, 3);
mov(1:nframe)=struct('cdata',[],'colormap',[]);
figure;
for k=1:nframe
    imagesc(data(:, :, k));
    mov(k)=getframe(gcf);
end
if nargin > 1 && ~isempty(filename)
    vidObj = VideoWriter(filename);
    open(vidObj);
    
    for k = 1:length(mov)
        
        writeVideo(vidObj,mov(k));
    end
    
    % Close the file.
    close(vidObj);
end