%% Play GCaMP movie, Thorlabs Body, Trajectory (PARALLELIZED!)

% make movie with this

% 2/17/16
% tic
clear all
cd('C:\Users\R123a\Desktop\ImagingAnalysisCode');
%%D:\GCAMPImaging\PostandPrettx3\GCaMP_03072022\03072022exp2.tif
global fps
fps = 15.625;
% fps = 30.51;
expDate = '03082022'; %date of experiment
fileFolders = {'D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp1\','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp2\','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp3\','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp4\'};
fileDirs={'D:\GCAMPImaging\PostandPrettx3\track_03082022\exp1\','D:\GCAMPImaging\PostandPrettx3\track_03082022\exp2\','D:\GCAMPImaging\PostandPrettx3\track_03082022\exp3\','D:\GCAMPImaging\PostandPrettx3\track_03082022\exp4\'};
Masks = {'D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp1mask.tif','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp2mask.tif','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp3mask.tif','D:\GCAMPImaging\PostandPrettx3\GCaMP_03082022\exp4mask.tif'};

for experiment=1:4
close all
gcp('nocreate');    %clear existing par pool
tic
expNum = ['exp',num2str(experiment)];   %experiment number
fileFolder = fileFolders{experiment};    %where you save GCaMP
fileDir=fileDirs{experiment}; %where you save track & track_FPGA

fileFolderTracker = fileDir;
addpath(fileFolder,fileFolderTracker);
toc
%load dir----------------------------------------
tic
testt=dir(fullfile(fileFolder,'*.tiff'));
clear dirOutput
firstFrame=length(testt);
for i = 1:length(testt)
    [b, remain]=strtok(testt(i).name, '_');
    b=strtok(remain, '_');
    b=str2double(b);
    dirOutput(b)=testt(i);
    if firstFrame>b
        firstFrame=b;
    end
end
clear b remain
fileNames = {dirOutput.name}';

i=1;
while isempty(fileNames{i}==1)
    i=i+1;
end
firstFrame=i;
fileNames=fileNames(firstFrame:end); 
    
toc
fprintf('done loading...\n');
%% Generate Tables (PARALLEL)
tic
poolobj = gcp; % If no pool, create new one. alternative: gcp('nocreate');
fprintf('generating tables...\n');
[FpgaTable, TrackTable] = getTable(fileDir); 
FpgaGCaMP=FpgaTable(:,3);   %FPGA corresp. to all GCaMP. units: GCaMPframe#
FpgaTracker=FpgaTable(:,2);
toc
fprintf('done generating tables...\n');
%% StartEndFrame,ExtractDatafromTables
tic
startFrame = max(200,FpgaGCaMP(1)+50); %GCaMP frame num of start (either frame 200 or 50th frame after saving)
% startFrame = 10000% 6885;
endFrame = firstFrame+length(fileNames)-200; %GCaMP frame num of end
%endFrame = 13000; %7350
nFrames = endFrame-startFrame+1;
frameOffset = startFrame-firstFrame+1;

GCaMPfrm = [startFrame:endFrame];   %GCaMP of interest

%find FPGA corresponding to GCaMP of interest
GCaMP2FPGA = zeros(1,nFrames);  %preallocate mem
GCaMP2Tracker = zeros(1,nFrames);   %preallocate mem
for i = 1:length(GCaMPfrm)
    tempFind = find(FpgaGCaMP>=GCaMPfrm(i),1);
    GCaMP2FPGA(i)=FpgaTable(tempFind,1);    %output FPGA frames corresp. to GCaMP of interest
    GCaMP2Tracker(i)=FpgaTable(tempFind,2);
    temp=tempFind;
end

% GCaMP2Tracker=FpgaTracker(GCaMP2FPGA);    %track frame corresp. to GCaMP of interest   
toc

% Extract data from Tables
[odor, CW, sync, time, servoX, servoY, piezoX, piezoY, bending,chr2,light]=getData(GCaMPfrm, FpgaTable, TrackTable);
trackX=servoX-3*piezoX; %um
trackY=servoY+3*piezoY; %um
[SmoothXY, Speed, AngSpeed, Reversals, ReversalsArray, Direction, Direction2,Time,trjX,trjY] = trackAnalysis(time, trackX, trackY);
SaveFileName = [expDate,'_',expNum,'_',num2str(startFrame),'-',num2str(endFrame),'_SmoothXY.mat'];
save(SaveFileName, 'SmoothXY');
SaveFileName2 = [expDate,'_',expNum,'_',num2str(startFrame),'-',num2str(endFrame),'_ReversalsArray.mat'];
save(SaveFileName2, 'ReversalsArray');
SaveFileName3 = [expDate,'_',expNum,'_',num2str(startFrame),'-',num2str(endFrame),'_odor.mat'];
save(SaveFileName3, 'odor');


%% Calc Max pixel value (USES PARALLEL POOL)
tic

frameDilution = 1;  %downsample n times
nFramesDil = floor(nFrames/frameDilution);
maxPixelDil = zeros(nFramesDil,1);
maxPixel = zeros(nFrames,1);
index = [];
maxPixelInd = zeros(nFrames,2);
maxPixelIndDil = zeros(nFramesDil,2);

mask = uint16(zeros(512,512));
    xLo = 271; yLo = 150;   %[X,Y]: values from imagesc(tempHeat');
    xHi = 400; yHi = 330;   %[X,Y]: values from imagesc(tempHeat');  
    mask(yLo:yHi, xLo:xHi) = 1;   %mask soma
if exist('BW2','var')
    BW = (BW2==0);
end
if exist('BW','var')
    mask = uint16(BW); %overwrite default mask
end

applyMask = 0;  %apply mask? true/false? !!!!!!!!!!!!!!!!!!!!!!!!!!

poolobj = gcp;
parfor t=1:nFramesDil
    frameNum = frameOffset+frameDilution*(t-1);
%     frameNum = frameOffset+(t-1)
    I2 = imread(fileNames{frameNum});
if applyMask==1 %apply mask True/False
    I2 = mask.*I2
end
    I2Old = I2;
%     imshow(2*I2Old);
%     pause(.2);
    I2 = reshape(I2,numel(I2),1);
    sorted = sort(I2(I2<62000),'descend');  %filter out cosmic rays
    index = [index, frameDilution*(t-1)+1]; %index corresp. to t
    maxPixelDil(t) = sorted(10);    %Take 10th brightest pixel
    [y, x] = find(I2Old==maxPixelDil(t),1); %pixel value may not be unique!
    maxPixelIndDil(t,:) = [x y]; 
    if mod(frameDilution*(t-1)+1,1000) <= frameDilution-1;
            fprintf('%d/%d frames maxed\n',frameDilution*(t-1)+1,nFrames); 
    end

end
delete(poolobj)

% map maxPixelDil to maxPixel using 'index'
for i=1:length(index)
    maxPixel(index(i))=maxPixelDil(i);
    maxPixelInd(index(i),:)=maxPixelIndDil(i,:);
end

for i=1:length(Reversals)   % get maxPixel at Reversals
    t = Reversals(i);
    frameNum = frameOffset+t-1;
    I2 = imread(fileNames{frameNum});
    I2Old = I2;
    I2 = reshape(I2,numel(I2),1);
    sorted = sort(I2(I2<62000),'descend');  %filter out cosmic rays
    maxPixel(t) = sorted(6);
    [y, x] = find(I2Old==maxPixelDil(t),1);
    maxPixelInd(t,:) = [x y]; 
end

toc
fprintf('done max pixel...\n');

%% save
SaveFileName = [expDate,'_',expNum,'_',num2str(startFrame), '-', num2str(endFrame), '_maxPixel.mat'];
save(SaveFileName, 'maxPixel','maxPixelInd');

SaveFileName = [expDate,'_',expNum,'_',num2str(startFrame), '-', num2str(endFrame), '_workspace.mat'];
save(SaveFileName);

fprintf('done saving workspace...\n');


BW2 = imread(Masks{1});
BW2 = (BW2==0);
figure; imshow(BW2)
%% Max Pixel 2
maxPixelDil2 = zeros(nFramesDil,1);
maxPixel2 = zeros(nFrames,1);
index2 = [];  
maxPixelInd2 = zeros(nFrames,2);
maxPixelIndDil2 = zeros(nFramesDil,2);

mask2 = uint16(zeros(512,512));
        xLo = 260; yLo = 180;   %[X,Y]: values from imagesc(tempHeat');
        xHi = 422; yHi = 370   ;   %[X,Y]: values from imagesc(tempHeat'); 
    mask2(yLo:yHi, xLo:xHi) = 1;   %mask soma
if exist('BW2','var')
    mask2 = uint16(BW2); %overwrite default mask
end
parfor t=1:nFramesDil
    frameNum = frameOffset+frameDilution*(t-1);
    I2 = imread(fileNames{frameNum});
% % %     xTemp = maxPixelInd(t,1); yTemp = maxPixelInd(t,2);
% % %     maskSize=35;
% % %     yLo = max(yTemp-maskSize,1); yHi = min(yTemp+maskSize,512); %set bounds
% % %     xLo = max(xTemp-maskSize,1); xHi = min(xTemp+maskSize,512); %set bounds
    I2 = mask2.*I2;
    I2Old = I2; %for finding index of max pixel
    I2 = reshape(I2,numel(I2),1); 
    sorted = sort(I2(I2<62000),'descend');  %filter out cosmic rays
    index2 = [index2, frameDilution*(t-1)+1]; %index corresp. to t
    maxPixelDil2(t) = sorted(6);
    [y, x] = find(I2Old==maxPixelDil2(t),1); %pixel value may not be unique!
    maxPixelIndDil2(t,:) = [x y]; 
%     imshow(I2Old*2)
%     hold on
%     plot(x,y,'r*');
%     pause
    if mod(frameDilution*(t-1)+1,1000) <= frameDilution-1;
            fprintf('%d/%d frames maxed\n',frameDilution*(t-1)+1,nFrames); 
    end
end

% map maxPixelDil to maxPixel using 'index'
for i=1:length(index)
    maxPixel2(index(i))=maxPixelDil2(i);
    maxPixelInd2(index(i),:)=maxPixelIndDil2(i,:);
end

fprintf('done max pixel2...\n');
figure;
MaxPixelAIY = maxPixel2
save(SaveFileName);
end
%% play movie
tic
fprintf('play movie. \n');
vidName = ['GCaMPBodyTrajMax_',expDate,'_',expNum,'.avi'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
% vid = VideoWriter('GCaMPBodyTrajMax_141110_exp5.avi','Motion JPEG AVI');
vid.Quality = 100;
vid.FrameRate = 24;
open(vid);

close(figure(1));
figure(1)
set(figure(1), 'Position', [10, 300, 720, 680]); %[20, 70, 1024, 620]
ymin2 = 4000;%min(maxPixel(maxPixel>0));
ymax2 = 20000;%max(maxPixel(maxPixel>0));
SmoothXYmm = SmoothXY/1000;

if exist('maxPeaks','var')
    tempValid=find(maxPeaks);
end
for t=1:20:nFrames-75%1:15:nFrames-75
    frameNum = frameOffset+t-1; %frame number according to the fileNames

subplot(3,3,1)  %plot GCaMP image
    hold off
    if exist('maxPeaks','var')
        [~,idx] = min(abs(tempValid-t));
        frameNum = frameOffset+tempValid(idx)-1;
    end
    I2 = imread(fileNames{frameNum}); %IMPORTANT! recall that data(t) corresponds to the best location we found at image t.  
    I2 = I2*2;%     I = imcrop(I,[600 550 400 400]);
    GCaMPframeNum = strsplit(fileNames{frameNum},'_');
    GCaMPframeNum = GCaMPframeNum{2};
    imshow(I2);
%     hold on
%     nearest = find(maxPeaks(1,t:end)>10000,1);
%     plot(maxPixelInd(t+nearest-1,1),maxPixelInd(t+nearest-1,2),'ro');
%     if size(maxPeaks,1)==2  %if E 2nd peak
%     nearest = find(maxPeaks(2,t:end)>10000,1);
%     plot(maxPixelInd(t+nearest-1,1),maxPixelInd(t+nearest-1,2),'bo');
%     end
%     nearest = find(maxPeaks2(1,t:end)>6000,1);
%     plot(maxPixelInd2(t+nearest-1,1),maxPixelInd2(t+nearest-1,2),'go');
    title(sprintf('GCaMP frm# %s',GCaMPframeNum),'FontWeight','normal','FontSize',12);
%      title('50X GCaMP','FontWeight','normal');
    hold off

subplot(3,3,4)  %plot Tracker image
    hold off
    trackerNum = GCaMP2Tracker(t);
%     bodyFileName = fileNamesTracker{find(fileNamesTrackerIndices==trackerNum,1),1};
%     Ibody = imread(fileNamesTracker{find(fileNamesTrackerIndices==trackerNum,1),1});
    tempBodyText = strcat('track_*_',num2str(trackerNum),'_*.jpg');
    tempBody=dir(fullfile(fileFolderTracker,tempBodyText));
    tempBody = tempBody.name;
    Ibody = imread(tempBody);
    Ibody = imcrop(Ibody,[500 380 640 640]);
    Ibody = imrotate(Ibody,90);
    imshow(Ibody);
    title(sprintf('tracker frm#: %d',trackerNum),'FontWeight','normal','FontSize',12);
%     title('1X Darkfield','FontWeight','normal');
    axis equal
    hold off

subplot(3,3,[2,3,5,6])  %plot Traj
    hold off
    if exist('lawnVector','var');
    pos = [SmoothXY(lawnEnter,1)-lawnVector(1)-lawnRadius SmoothXY(lawnEnter,2)+lawnVector(2)-lawnRadius 2*lawnRadius 2*lawnRadius]/1000;
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[.95 .95 .95]);
    hold on
    end
    plot(SmoothXYmm(1:10:t,1),SmoothXYmm(1:10:t,2), 'Color','blue');  %downsample plot by 10 for speed
    hold on
    plot(SmoothXYmm(1,1),SmoothXYmm(1,2),'kx');
    plot(SmoothXYmm(t:10:end,1),SmoothXYmm(t:10:end,2), 'Color',[.6 .6 .6]);
% % % %     plot(SmoothXYmm(Reversals,1),SmoothXYmm(Reversals,2),'Color',[.6 .6 .6],'Marker','.','LineStyle','none');
% % % %     plot(SmoothXYmm(Reversals(Reversals<=t),1),SmoothXYmm(Reversals(Reversals<=t),2),'m.');
    hold off
    axis([SmoothXYmm(t,1)-7 SmoothXYmm(t,1)+7 SmoothXYmm(t,2)-7 SmoothXYmm(t,2)+7]); 
    title('Trajectory','FontWeight','normal');
    
subplot(3,3,[7,8,9])  %plot GCaMP Intensity
    hold off
    if t==1
        if exist('touchON','var')
            plot((startFrame+(1:nFrames))/(fps*60),touchON*ymax2,'Color',[.7 .7 .7])
            hold on
        end
        if max(odor>0)
            area((startFrame:endFrame)/(fps*60),ymax2*odor,'FaceColor',[.2 .2 .2],...
            'EdgeColor',[.2 .2 .2]);
            hold on
        end
        if exist('maxPeaksInterp','var')
            xVal = (startFrame+(1:length(maxPeaksInterp)))/(fps*60);
            yVal = maxPeaksInterp(1,:);
        else
            xVal =(startFrame+index-1)/(fps*60);    %index=[1,1+frameDil,1+2*frameDil,...]
            yVal = maxPixelDil;
        end
        xVal = xVal(1:30:end); yVal = yVal(1:30:end);   %downsample for speed
        plot(xVal,yVal,'b')
        hold on
        if size(maxPeaksInterp,1)>1
            yVal = maxPeaksInterp(2,:); yVal = yVal(1:30:end);
            plot(xVal,yVal,'r')
        end
        hLine = line([(startFrame+t-1)/(fps*60) (startFrame+t-1)/(fps*60)],[ymin2 ymax2],'Color','green');
        xlabel('time (min)'); ylabel('max pixel intensity');
        axis([(startFrame)/(fps*60) (endFrame)/(fps*60) ymin2 ymax2]);
    end
    hLine.XData = [(startFrame+t-1)/(fps*60) (startFrame+t-1)/(fps*60)];
%     plot((startFrame+Reversals)/(fps*60),ymax2/3,'m.','MarkerSize',4);
    hold off
    writeVideo(vid, getframe(figure(1)));   %getframe(gcf);
end

close(vid);
fprintf('done making movie\n');
toc

%% Make figure (Trajectory w GCaMP color, AIY timeseries w odor & reversals)
h=figure(2);
% set(figure(1), 'Position', [10, 300, 2500, 640]); %[20, 70, 1024, 620]
set(figure(2),'Position',[10,10,768,1080]);
% set(h, 'units','normalized','outerposition', [0 0 1 1])
ymin2 = 0%6000;%min(maxPixel(maxPixel>0));
ymax2 = 255%40000;%max(maxPixel(maxPixel>0));

%trajectory w GCaMP color
% subplot(1,4,[1])
subplot(3,2,[1 2 3 4]);
hold on
colormap(jet)

if exist('lawnVector','var');
    lawnEnter = find(abs(GCaMP2Tracker-trackerNum)<=1,1);   %units frame
pos = [SmoothXY(lawnEnter,1)-lawnVector(1)-lawnRadius SmoothXY(lawnEnter,2)+lawnVector(2)-lawnRadius 2*lawnRadius 2*lawnRadius]/1000;
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[.95 .95 .95]);
end

if exist('maxPeaksInterp','var')
    scatter(SmoothXY(:,1)/1000,SmoothXY(:,2)/1000,5,maxPeaksInterp(1,:),'.')
else
    maxPixelInterp = interp1(index,maxPixelDil,startFrame:endFrame);    %wrong!
    maxPixelInterp = interp1(index,maxPixelDil,index(1):index(end));
    GCaMPThresh = .8*max(maxPixelInterp);   %80% of max
    maxPixelInterp(maxPixelInterp>GCaMPThresh)=GCaMPThresh;
    scatter(SmoothXY(:,1)/1000,SmoothXY(:,2)/1000,5,maxPixelInterp(1,:),'.')
end
caxis([0 210]) %
plot(SmoothXY(1,1)/1000,SmoothXY(1,2)/1000,'ko');   %'o' @ start
plot(SmoothXY(end,1)/1000,SmoothXY(end,2)/1000,'kx');   %'x' @ end
%add time marks
% for i=1:940:length(SmoothXY(:,1))
for i=1:500:length(SmoothXY(:,1))
    text(SmoothXY(i,1)/1000,SmoothXY(i,2)/1000,num2str(round(10*(startFrame+i)/(fps*60))/10))
end
axis equal
xlabel('mm'); ylabel('mm');
title(['Trajectory']);
colorbar;

hold off

%AIY timeseries w odor % reversals
% subplot(1,4,[2 3 4])
subplot(3,2,[5 6])
area((startFrame:endFrame)/(fps*60),ymax2*odor,'FaceColor',[.3 .3 .3],'EdgeColor',[.3 .3 .3]);
hold on
if exist('touchON','var')
    plot((startFrame+(1:nFrames))/(fps*60),touchON*ymax2,'Color',[.6 .6 .6])
end
if exist('maxPeaksInterp','var')
    xVal = (startFrame+(1:length(maxPeaksInterp)))/(fps*60);
    yVal = maxPeaksInterp(1,:);
else
    xVal =(startFrame+index-1)/(fps*60);    %index=[1,1+frameDil,1+2*frameDil,...]
    yVal = maxPixelDil;
end
plot(xVal,yVal,'b')
if exist('maxPeaksInterp','var')
    if size(maxPeaksInterp,1)>1
        yVal = maxPeaksInterp(2,:); yVal = yVal(1:30:end);
        plot(xVal(1:30:end),yVal,'r')
    end
end
% plot((startFrame+Reversals)/(fps*60),25000,'m.','MarkerSize',4);
% axis([(startFrame)/(fps*60) (endFrame)/(fps*60) ymin2 ymax2]);
% axis([xVal(1) 1.02*xVal(end) ymin2 ymax2]);
% axis([(startFrame)/(fps*60) 10 ymin2 ymax2]);    %PLOT UP TO 60mins
xlabel('time (min)'); ylabel('max pixel intensity');
title(['GCaMP activity ',expDate,'_.',expNum]);
hold off
figName = [expDate,'_',expNum,'_GCaMPTraj.fig'];
figNamePNG = [expDate,'_',expNum,'_GCaMPTraj.png'];
savefig(h,figName)
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r0',figNamePNG)
fprintf('done making figure \n');    


end

%% find max per stack
nSlices = 30;
maxSlice = [];
test=[];
maxPosition = [];
phi = 1;    %offset z-scan
for i=phi:nSlices:length(maxPixel)-nSlices
    [maxTemp,maxTempIndex] = max(maxPixel(i:i+(nSlices-1)));
    maxSlice(i+maxTempIndex-1)=maxTemp;
    test(i)=1;
    maxPosition(i) = maxTempIndex;
end
figure;
plot(maxPixel);
hold on
plot(maxSlice)
plot(test*max(maxSlice))

%% find n greatest local maxima per stack
locArray = [];
nLocalMaxToFind = 1;
for i=phi:nSlices:length(maxPixel)-nSlices
    [peak,loc] = findpeaks(maxPixel(i:i+(nSlices-1)));
    if length(peak)>=nLocalMaxToFind
        peakSort = sort(peak,'descend');
        for i=1:nLocalMaxToFind
            tempFind = find(peak==peakSort(i));
            locArray = [locArray; loc(tempFind)];
        end
    end
end
figure
histogram(locArray,[0.5:1:nSlices]); title(['hist of ',...
    num2str(nLocalMaxToFind),' largest local maxima'])
%%
ranges = [];
maxPeaks = [];
max_1 = []; max_2 = [];
goodpeak = [];
searchRange = [1 30]; %[lo1 hi1; lo2 hi2; ...]
% % % % % % % lo_1 = 5; hi_1 = 11;
% % % % % % % lo_2 = 11; hi_2 = 18;
for i=phi:nSlices:length(maxPixel)-nSlices
    for j=1:nLocalMaxToFind
        lo = searchRange(j,1); hi = searchRange(j,2);
        [peak, loc] = findpeaks(maxPixel(i+lo-1:i+hi-1));
        if isempty(loc) %set=0 if you cannot find peak (poor resolution);use local max if you cannot find peak
%         [maxTemp,maxTempIndex] = max(maxPixel(i+lo_1-1:i+hi_1-1));
%         max_1(i+lo_1+maxTempIndex-1-1)=maxTemp;
        else
            [maxPk, maxPkInd] = max(peak);
            maxPeaks(j,i+lo+loc(maxPkInd)-1-1)=peak(maxPkInd);%peaks(maxPkInd);
            goodpeak(i+lo+loc(maxPkInd)-1-1)=peak(maxPkInd);%peaks(maxPkInd);
        end
    end
%     ranges(i+lo_2-1)=20000;
%     ranges(i+hi_2-1)=20000;
end
%pad maxPeaks to be length of maxPixel
if length(maxPeaks)<length(maxPixel)
    maxPeaks(end:length(maxPixel))=0;
    maxPeaks(end)=maxPixel(end);
end
maxPeaks(1)=maxPixel(1);
tempIndex = find(maxPeaks);
maxPeaksInterp = interp1(tempIndex,maxPeaks(tempIndex),1:length(maxPeaks));

figure; 
plot(maxPixel);
hold on
plot(maxPeaks(1,:),'*');
% plot(maxPeaks(2,:),'o');
plot(test*.9*max(maxSlice));

%% FFT to get period
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
data = maxPixel;
data = data(1:min(4096,length(data)));
L = length(data);
t = (0:L-1)*T;                % Time vector
y=data;
% y=data.*hamming(L)'; %apply Hamming window
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% NFFT = L;   %don't care about NLogN vs N^2 speed
Y = fft(y,NFFT)/L;
% % f = Fs/2*linspace(0,1,NFFT/2+1);
f = 1/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
amp=2*abs(Y(1:NFFT/2+1));
figure; plot(f,10*log(2*abs(Y(1:NFFT/2+1)))) 
% figure; plot(f,2*abs(Y(1:NFFT/2+1)));
title('Single-Sided Amplitude Spectrum of y(t) dB')
xlabel('Frequency (Hz)')
ylabel('|Y(f)| dB')
% temp=find(f>.033,1); %30 slices per stack max
% amp(1:temp)=0;
% tempIndex = find(amp==max(amp));
% fprintf('test');
% nSlices=round(1/f(tempIndex));

%% Detect Lawn

fileFolderTracker = fileDir;
% trackerNum = GCaMP2Tracker(GCaMPLawn);
    trackerNum = 34218; %tracker frame showing lawn (prob second #)
    tempBodyText = strcat('track_*_',num2str(trackerNum),'_*.jpg');
    tempBody=dir(fullfile(fileFolderTracker,tempBodyText));
    tempBody = tempBody.name
    Ibody = imread([fileFolderTracker,tempBody]);
%     Ibody = imcrop(Ibody,[360 240 920 920]);
    Ibody = imrotate(Ibody,90);
    h=figure(10);
    imshow(Ibody);
    title('click ENTRY PT FIRST then TWO other pts');
    %title(sprintf('tracker frame#: %d',trackerNum),'FontWeight','normal');
%     title('1X Darkfield','FontWeight','normal');
    [xs,ys] = ginput(3);
    [xfit,yfit,Rfit] = circfit(xs,ys);
    hold on
    rectangle('Position',[xfit-Rfit,yfit-Rfit,2*Rfit,2*Rfit],'Curvature',[1 1]);
    pause(1)
    close(h)
    lawnVector = 5.168*[xs(1)-xfit, ys(1)-yfit];
    R = [cosd(1) -sind(1); sind(1) cosd(1)]; lawnVector=lawnVector*R;   %account for 3deg angle between thor camera and stage
    lawnRadius = 5.168*Rfit;

%% Trim Trajectory 
if ~exist('GCaMPLawn','var')
    GCaMPLawn = 1;
end
[trimStart,trimEnd] = trimTrajectory(SmoothXY,Speed,GCaMPLawn);
SmoothXY2=SmoothXY(trimStart:trimEnd,:); Speed2=Speed(trimStart:trimEnd); GCaMPLawn2 = GCaMPLawn-trimStart;
    
%% Touch detection
tic
[touchON,testOut] = detectTouch(GCaMP2Tracker,fileFolderTracker,3);
toc

%% Widefield movie for Hannah
tic
fprintf('play movie. \n');
vidName = ['Body_',expDate,'_',expNum,'.avi'];
vid = VideoWriter(vidName,'Motion JPEG AVI');
% vid = VideoWriter('GCaMPBodyTrajMax_141110_exp5.avi','Motion JPEG AVI');
vid.Quality = 100;
vid.FrameRate = 24;
open(vid);

figure(1)

for t=1:2:nFrames
    trackerNum = GCaMP2Tracker(t);
%     bodyFileName = fileNamesTracker{find(fileNamesTrackerIndices==trackerNum,1),1};
%     Ibody = imread(fileNamesTracker{find(fileNamesTrackerIndices==trackerNum,1),1});
    tempBodyText = strcat('track_*_',num2str(trackerNum),'_*.jpg');
    tempBody=dir(fullfile(fileFolderTracker,tempBodyText));
    tempBody = tempBody.name;
    Ibody = imread(tempBody);
    Ibody = imcrop(Ibody,[620 500 400 400]);
    Ibody = imrotate(Ibody,90);
    imshow(Ibody);
    %title(sprintf('tracker frame#: %d',trackerNum),'FontWeight','normal');
    title('Gradual Turning','FontWeight','normal');    
    writeVideo(vid, getframe(figure(1)));   %getframe(gcf);

%    pause(.05);
end
% 
close(vid);
fprintf('done making movie\n');
toc

%% INDEX of VARIABLES