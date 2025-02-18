% BME345_Project_pt5.m
% 02/10/2025
%
% Track and make stick figure of cycling motion

clear
clc
close all


%% Declarations
% Video Info
vidFile = 'Bike_Mirrored.mov';
frameStart = 1; 
frameStop = 60; 
vid = VideoReader(vidFile);

% Figure Background
background = 'Bike_Background.png';

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

r1 = 0.665; %hip to bike pedal (m)
r2 = 0.185; %bike pedal to foot (m)
r3 = 0.44; %foot to knee (m)
r4 = 0.51; %knee to hip (m)

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

%% Pixel to Meter Conversion & Vector Math for Angle Guesses

%Pixel position vectors
r1Vec = [YxPos(1)-GxPos(1), YyPos(1)-GyPos(1)];
r2Vec = [BxPos(1)-YxPos(1), ByPos(1)-YyPos(1)];
r3Vec = [RxPos(1)-BxPos(1), RyPos(1)-ByPos(1)];
r4Vec = [GxPos(1)-RxPos(1), GyPos(1)-RyPos(1)];

%Pixel magnitudes
r1PxMag = vecnorm(r1Vec);
r2PxMag = vecnorm(r2Vec);
r3PxMag = vecnorm(r3Vec);
r4PxMag = vecnorm(r3Vec);

pxmconv = [r1/r1PxMag, r2/r2PxMag, r3/r3PxMag, r4/r4PxMag];
px_mConv = mean(pxmconv); % Conversion factor for px to m

% Convert position vectors to m
RxPos = RxPos.*px_mConv;
RyPos = RyPos.*px_mConv;

GxPos = GxPos.*px_mConv;
GyPos = GyPos.*px_mConv;

BxPos = BxPos.*px_mConv;
ByPos = ByPos.*px_mConv;

YxPos = YxPos.*px_mConv;
YyPos = YyPos.*px_mConv;

% figure; %Four-bar position at very start of video
% plot(GxPos(1), GyPos(1), 'g.', ...
%                 RxPos(1), RyPos(1), 'r.', ...
%                 BxPos(1), ByPos(1), 'c.', ...
%                 YxPos(1), YyPos(1), 'y.', ...
%                 [GxPos(1), YxPos(1)], [GyPos(1), YyPos(1)], ...
%                 [YxPos(1), BxPos(1)], [YyPos(1), ByPos(1)],  ...
%                 [BxPos(1), RxPos(1)], [ByPos(1), RyPos(1)], ...
%                 [RxPos(1), GxPos(1)], [RyPos(1), GyPos(1)], ...
%                 'LineWidth', 1, ...
%                 'MarkerSize', 20)
% title("Starting Position")

% Convert videoframe1 r vectors units to m
r1Vec = r1Vec.*px_mConv;
r2Vec = r2Vec.*px_mConv;
r3Vec = r3Vec.*px_mConv;
r4Vec = r4Vec.*px_mConv;

% +x axis to find videoframe1 angles 
xAxisVec = [1, 0];
xAxisMag = vecnorm(xAxisVec);

% Angles in radians
th1 = (2*pi) - acos(dot(r1Vec, xAxisVec) ./ (r1 .* xAxisMag)); %adjusted so angle is in CCW direction
th2 = acos(dot(r2Vec, xAxisVec) ./ (r2 .* xAxisMag));
th3 = acos(dot(r3Vec, xAxisVec) ./ (r3 .* xAxisMag));
th4 = acos(dot(r4Vec, xAxisVec) ./ (r4 .* xAxisMag)); 


%% Q5: Stick Figure Video

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
        xlim([0, col*px_mConv])
        ylim([0, row*px_mConv])
        pbaspect([1, row/col, 1])  

    drawnow % forces figure to appear, which may not happen in loops 

    frame = getframe(gcf);
    writeVideo(vidObj2,frame); % writes frame to video

end

close(vidObj2) % closes video
