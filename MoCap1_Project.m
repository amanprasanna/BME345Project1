% Mylah Williams
% mlwill28@ncsu.edu
% 02/10/2025
% MoCap1_Project.m
%
% Tracks LEDs' motion and makes an animated stick figure of the motion

clear
clc
close all


%% Declarations
% Video Info
vidFile = 'IMG_0217.mov';
frameStart = 1; 
frameStop = 185; 
vid = VideoReader(vidFile);

% For cropping
widthEdge = 700;
widthRange = 600;
heightEdge = 1;
heightRange = 840;

% Color Adjustments
lthrR = 220; %red tape lower threshold
uthrR = 255; %red tape upper threshold

lthrG = 220; %green tape lower threshold
uthrG = 255; %green tape upper threshold

lthrRB = 106;%red lthr for blue tape
uthrRB = 120;%red uthr for blue tape
lthrGB = 155;%green lthr for blue tape
uthrGB = 170;%green uthr for blue tape
lthrBB = 185;%blue lthr for blue tape
uthrBB = 210;%blue uthr for blue tape

lthrGen = 110; %generalized lower threshold
uthrGen = 225; %generalized upper threshold

lthrRY = 220;%red lthr for yellow tape
lthrGY = 170;%green lthr for yellow tape
lthrBY = 40;%blue lthr for yellow tape
uthrBY = 85;%blue uthr for yellow tape

%% Video 1: Centroid Movement

vidObj1 = VideoWriter("MoCapP1V1_Williams");
vidObj1.FrameRate = vid.FrameRate; % set frame rate (e.g., 15)
open(vidObj1) % opens video for recording 

% Step through each frame 
for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice
    frameSlice = flip(frameSlice, 2);

    frameSliceCrop = imcrop(frameSlice, ...
        [widthEdge, heightEdge, widthRange, heightRange]);
    red = frameSliceCrop(:,:,1);
    green = frameSliceCrop(:,:,2);
    blue = frameSliceCrop(:,:,3);
    RframeSliceBW = ((red >= lthrR) & (red <= uthrR) & ...
                    (green <= lthrGen) & (blue <= lthrGen)); 
    GframeSliceBW = ((green >= lthrG) & (green <= uthrG) & ...
                    (red <= lthrGen) & (blue <= lthrGen)); 
    BframeSliceBW = ((blue >= lthrBB) & (blue <= uthrBB)) & ...
                    ((red >= lthrRB) & (red <= uthrRB) & ...
                    (green >= lthrGB) & (green <= uthrGB)); 
    YframeSliceBW = ((red >= lthrRY) & (red <= uthrR) & ...
                    (green >= lthrGY) & (green <= uthrG) & ...
                    (blue >= lthrBY) & (blue <= uthrBY));

    frameSliceBW = RframeSliceBW + GframeSliceBW + BframeSliceBW + YframeSliceBW;
    [row,col, ~] = size(frameSliceBW);

    % Uses custom function centroid to return row and col of a binary image
    [RcentRow(k), RcentCol(k)] = Centroid345(RframeSliceBW); %#ok<SAGROW>
    [GcentRow(k), GcentCol(k)] = Centroid345(GframeSliceBW); %#ok<SAGROW>
    [BcentRow(k), BcentCol(k)] = Centroid345(BframeSliceBW); %#ok<SAGROW>
    [YcentRow(k), YcentCol(k)] = Centroid345(YframeSliceBW); %#ok<SAGROW>
    
    % Returns x and y positions of centroids
    RxPos = RcentCol(frameStart:end);
    RyPos = row - RcentRow(frameStart:end);

    GxPos = GcentCol(frameStart:end);
    GyPos = row - GcentRow(frameStart:end);

    BxPos = BcentCol(frameStart:end);
    ByPos = row - BcentRow(frameStart:end);

    YxPos = YcentCol(frameStart:end);
    YyPos = row - YcentRow(frameStart:end);

    % Images + Plot Output
    figure (1)
        subplot(3, 1, 1)
        imshow(frameSliceCrop)
        title('Color Cropped Image')
        pbaspect([1, row/col, 1])

    subplot(3, 1, 2)
        imshow(frameSliceBW)
        title('Thresholded Image')
        pbaspect([1, row/col, 1])

    subplot(3, 1, 3)
        for h = 1:length(GxPos)
          plot(RxPos, RyPos, 'r', RxPos(h), RyPos(h), 'rx-', ...
              GxPos, GyPos, 'g', GxPos(h), GyPos(h), 'gx-', ...
              BxPos, ByPos, 'c', BxPos(h), ByPos(h), 'cx-', ...
              YxPos, YyPos, 'y', YxPos(h), YyPos(h), 'yx-', 'LineWidth', 1)
        end
        title('Joint Centroids')
        xlim([0, col])
        ylim([0, row])
        set(gca, 'XTick', [], 'YTick', [])
        pbaspect([1, row/col, 1])
    drawnow % forces figure to appear, which may not happen in loops 

    frame = getframe(gcf);
    writeVideo(vidObj1,frame); % writes frame to video 

end

close(vidObj1) % closes video


%% Video 2: Stick Figure

t = (length(frameStart:frameStop))/vid.FrameRate;
    tLin = linspace(0, t, length(frameStart:frameStop));

    % Velocity in px/s
    vRx = gradient(RxPos, tLin);
    vRy = gradient(RyPos, tLin);

    vGx = gradient(GxPos, tLin);
    vGy = gradient(GyPos, tLin);

    vBx = gradient(BxPos, tLin);
    vBy = gradient(ByPos, tLin);

    vYx = gradient(YxPos, tLin);
    vYy = gradient(YyPos, tLin);

vidObj2 = VideoWriter("MoCapP1V2_Williams");
vidObj2.FrameRate = vid.FrameRate; % set frame rate (e.g., 15)
open(vidObj2) % opens video for recording 

for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice

    frameSliceCrop = imcrop(frameSlice, ...
        [widthEdge, heightEdge, widthRange, heightRange]);

    % Images + Plot Output
    figure (2)
        subplot(2, 1, 1)
        imshow(frameSliceCrop)
        title('Color Cropped Image')
        pbaspect([1, row/col, 1])

    subplot(2, 1, 2)
        for h = 1:(k-frameStart)
            plot(RxPos(h), RyPos(h), 'r.', ...
                GxPos(h), GyPos(h), 'g.', ...
                BxPos(h), ByPos(h), 'c.', ...
                YxPos(h), YyPos(h), 'y.', ...
                [RxPos(h), GxPos(h)], [RyPos(h), GyPos(h)], 'r-', ...
                [GxPos(h), BxPos(h)], [GyPos(h), ByPos(h)], 'g-', ...
                [BxPos(h), YxPos(h)], [ByPos(h), YyPos(h)], 'c-', ...
                'LineWidth', 1, ...
                'MarkerSize', 25)
            hold on
            quiver(RxPos(h), RyPos(h), vRx(h), vRy(h), .1, "Color", "r", ...
                "LineWidth", 2, "MaxHeadSize", 1)
            quiver(GxPos(h), GyPos(h), vGx(h), vGy(h), .1, "Color", "g", ...
                "LineWidth", 2, "MaxHeadSize", 1)
            quiver(BxPos(h), ByPos(h), vBx(h), vBy(h), .1, "Color", "c", ...
                "LineWidth", 2, "MaxHeadSize", 1)
            quiver(YxPos(h), YyPos(h), vYx(h), vYy(h), .1, "Color", "y", ...
                "LineWidth", 2, "MaxHeadSize", 1)
            hold off
        end
        title('Joint Velocities (px/s)')
        xlim([0, col])
        ylim([0, row])
        set(gca, 'XTick', [], 'YTick', [])
        pbaspect([1, row/col, 1])
    drawnow % forces figure to appear, which may not happen in loops 

    frame = getframe(gcf);
    writeVideo(vidObj2,frame); % writes frame to video

end

close(vidObj2) % closes video
