deltaFOverFData=readmatrix('/misc/public/For David/∆FF vs time.csv');
timeS = deltaFOverFData(:,1);
deltaFOverF = 100*deltaFOverFData(:,2:end)';

%V = tiffreadVolume('/misc/public/For David/ROI-Composite.tif');
fireLUT = readmatrix('fireLUT.csv');
fireLUT = fireLUT(:,2:end)/255.0;

numColors = 256;
redColormap = [linspace(0, 1, numColors)', zeros(numColors, 1), zeros(numColors, 1)];
% red: 0,16464
% df: 100, 600
%colormap(fireLUT(:,2:end)/255.0)
deltaFOverFmean = mean(deltaFOverF);
deltaFOverFstd = std(deltaFOverF);
for timestep=1:length(timeS)
    close all;
    [neuronImg,deltaFOverFim]=combineImages(V(:,:,2*(timestep-1)+2), V(:,:,2*(timestep-1)+1));
    figure('Color', 'k',"Units","inches","Position",[0,0,4.5,2.47],'Resize','off','Visible','off');
    tiledlayout(1,2,"TileSpacing","tight","Padding","compact")%,"Units","inches","OuterPosition",[0,0,4.63,2.47]);
    set(gcf, 'InvertHardcopy', 'off');
    set(gcf, 'PaperPositionMode', 'auto');

    nexttile
    imshow(neuronImg+deltaFOverFim)
    
    %% 20 um scale bar
    % for now hardcode pixel size
    pixelSize = 0.5405;

    axis equal
    axis square
    rectangle('Position',[105,138,20/pixelSize,3],'LineWidth',2,'LineStyle','none','FaceColor','w')
    % time
    txt = sprintf('00:%02.f sec',floor(timeS(timestep)));
    text(5,15,txt,"color","w",'FontSize',14,'FontWeight','bold')

    %subplot(1,2,2);

    
    % create a figure with two subplots side by side and with a black background
    nexttile
    set(gca,'Color','k')
    plotColor = [197/255.0, 66/255.0, 245/255.0];
    if timestep<length(timeS)
        plotWithShadedError(timeS(timestep+1:end),deltaFOverFmean(timestep+1:end), deltaFOverFstd(timestep+1:end), 'w');
    end
    plotWithShadedError(timeS(1:min(timestep+1,length(timeS))),deltaFOverFmean(1:min(timestep+1,length(timeS))), deltaFOverFstd(1:min(timestep+1,length(timeS))), plotColor);
    ylim([-.01,.21]*100)
    xlim([0 55])
    set(gca,'XColor', 'w','YColor','w')
    ylabel('{\DeltaF}/F (%)')
    xlabel("Time (s)")
    %text(20,.18,"\DeltaF/F","color",plotColor,'FontSize',14,'FontWeight','bold')
    axis square

    % rescale to keep the same size
    %h1 = subplot(1, 2, 1); % Handle for the first subplot
    %h2 = subplot(1, 2, 2); % Handle for the second subplot
    
    % Set the positions to be the same
    % pos1 = get(h1, 'Position');
    % pos2 = get(h2, 'Position');
    
    % Use the smaller width and height for both subplots
    newWidth = .3347;%max(pos1(3), pos2(3));
    newHeight = .3347;%max(pos1(4), pos2(4));
    
    % Set the new positions
    %set(h1, 'Position', [pos1(1), pos1(2), newWidth, newWidth]);
    %set(h2, 'Position', [pos2(1)-.15, pos2(2), newWidth, newWidth]);
    %set(gcf,'Color',[0 0 0]);
    outputName = sprintf("/groups/beyene/beyenelab/Imaging Data/Ackerman/compositeMovie/tifs/%03d.tif",timestep);
    %exportgraphics(gcf, outputName, 'BackgroundColor', 'k','Resolution',300);
    print(gcf, outputName, '-dtiff', '-r300');
    %saveas(gcf,'temp.tif','Resolution',300);
end



% after running the above, compile it into mov: ffmpeg -r 30 -f image2 -i %03d.tif -vf "crop=1300:740:100:0,fps=30" -c:v libx264 -pix_fmt yuv420p -r 30 composite.mov