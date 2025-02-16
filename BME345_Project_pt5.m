% BME345_Project_pt5.m
% 02/10/2025
%
% Track and make stick figure of cycling motion

clear
clc
close all


%% Declarations
% Video Info
vidFile = 'IMG_0217.mov';
frameStart = 1; 
frameStop = 60; 
vid = VideoReader(vidFile);

% For cropping
widthEdge = 700;
widthRange = 630;
heightEdge = 10;
heightRange = 830;

% Color Adjustments
lthrR = 220; %red tape lower threshold
uthrR = 255; %red tape upper threshold

lthrRG = 135;%red lthr for green tape
uthrRG = 170;%red uthr for green tape
lthrBG = 140;%blue lthr for green tape
uthrBG = 190;%blue uthr for green tape
lthrGG = 220; %green tape lower threshold
uthrGG = 255; %green tape upper threshold

lthrRB = 106;%red lthr for blue tape
uthrRB = 120;%red uthr for blue tape
lthrGB = 155;%green lthr for blue tape
uthrGB = 170;%green uthr for blue tape
lthrBB = 185;%blue lthr for blue tape
uthrBB = 210;%blue uthr for blue tape

lthrGen = 110; %generalized lower threshold
uthrGen = 225; %generalized upper threshold

lthrRY = 185;%red lthr for yellow tape
uthrRY = 210;%red uthr for yellow tape
lthrGY = 170;%green lthr for yellow tape
uthrGY = 190;%green uthr for yellow tape
lthrBY = 0;%blue lthr for yellow tape
uthrBY = 10;%blue uthr for yellow tape

%% Determining Centroid Movement


% Step through each frame 
for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice

    frameSliceCrop = imcrop(frameSlice, ...
        [widthEdge, heightEdge, widthRange, heightRange]);
    red = frameSliceCrop(:,:,1);
    green = frameSliceCrop(:,:,2);
    blue = frameSliceCrop(:,:,3);
    RframeSliceBW = ((red >= lthrR) & (red <= uthrR) & ...
                    (green <= lthrGen) & (blue <= lthrGen)); 
    GframeSliceBW = ((green >= lthrGG) & (green <= uthrGG)) & ...
                    ((red >= lthrRG) & (red <= uthrRG)) & ...
                    ((blue >= lthrBG) & (blue <= uthrBG));
    BframeSliceBW = ((blue >= lthrBB) & (blue <= uthrBB)) & ...
                    ((red >= lthrRB) & (red <= uthrRB) & ...
                    (green >= lthrGB) & (green <= uthrGB)); 
    YframeSliceBW = ((red >= lthrRY) & (red <= uthrRY) & ...
                    (green >= lthrGY) & (green <= uthrGY) & ...
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

end

%% Stick Figure Video

t = (length(frameStart:frameStop))/vid.FrameRate;
    tLin = linspace(0, t, length(frameStart:frameStop));

% Force at each joint in N
% fGx = ;
% fGy = ;
% 
% fRx = ;
% fRy = ;
% 
% fBx = ;
% fBy = ;
% 
% fYx = ;
% fYy = ;

vidObj2 = VideoWriter("cycling");
vidObj2.FrameRate = vid.FrameRate; % set frame rate (e.g., 15)
open(vidObj2) % opens video for recording 

for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice

    frameSliceCrop = imcrop(frameSlice, ...
        [widthEdge, heightEdge, widthRange, heightRange]);

    % Images + Plot Output
    figure (2)
        for h = 1:(k-frameStart)
            plot(GxPos(1), GyPos(1), 'g.', ...
                RxPos(h), RyPos(h), 'r.', ...
                BxPos(h), ByPos(h), 'c.', ...
                YxPos(1), YyPos(1), 'y.', ...
                [GxPos(1), YxPos(1)], [GyPos(1), YyPos(1)], ...
                [YxPos(1), BxPos(h)], [YyPos(1), ByPos(h)],  ...
                [BxPos(h), RxPos(h)], [ByPos(h), RyPos(h)], ...
                [RxPos(h), GxPos(1)], [RyPos(h), GyPos(1)], ...
                RxPos, RyPos, 'r-', ...
                'LineWidth', 1, ...
                'MarkerSize', 20)             % CAN ADD CoM POINTS

            
            % CAN ADD FORCE ARROWS HERE
            % hold on
            % quiver(RxPos(h), RyPos(h), fRx(h), fRy(h), .1, "Color", "r", ...
            %     "LineWidth", 2, "MaxHeadSize", 1)
            % quiver(GxPos(h), GyPos(h), fGx(h), fGy(h), .1, "Color", "g", ...
            %     "LineWidth", 2, "MaxHeadSize", 1)
            % quiver(BxPos(h), ByPos(h), fBx(h), fBy(h), .1, "Color", "c", ...
            %     "LineWidth", 2, "MaxHeadSize", 1)
            % quiver(YxPos(h), YyPos(h), fYx(h), fYy(h), .1, "Color", "y", ...
            %     "LineWidth", 2, "MaxHeadSize", 1)
            % hold off
        end
        title('Pedaling System')
        xlabel('Horizontal Position (m)')
        ylabel('Vertical Position (m)')

        % CONVERT AXIS DIMENSIONS FROM PX TO METERS 

        legend('Hip', 'Knee', 'Foot', 'Bottom Bracket', 'Frame (r1)', ...
            'Pedal (r2)', 'Leg (r3)', 'Thigh (r4)', 'Location','eastoutside')
        xlim([0, col])
        ylim([0, row])
        set(gca, 'XTick', [], 'YTick', [])
        pbaspect([1, row/col, 1])
        % ADD STATIC BACKGROUND IMAGE
        % hold on 
        % I = imread('____.png');
        % h = image(xlim, -ylim, I);
        % uistack(h, 'bottom');

    drawnow % forces figure to appear, which may not happen in loops 

    frame = getframe(gcf);
    writeVideo(vidObj2,frame); % writes frame to video

end

close(vidObj2) % closes video