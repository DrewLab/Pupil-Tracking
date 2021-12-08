function[PupilTracker,BW,intensityThresh]=PupilTracker_QZ_001(cameraFs,nFramesToRead,pupilCamFileID,BW,intensityThresh,varargin)
%Inputs:
%pupilCamFileID: character string of .bin file of eye camera
%intensityThresh: pixel intensity threshold for pre-radon image
%thresholding to remove sclera
%BW: existing ROI for eye
%Outputs:
%PupilTracker.Pupil_Area: framewise area measurement in pixels
%PupilTracker.Overlay: movie of raw image data with pupil tracking data
%overlayed as blue mask.
%BW: Pre-defined region of the open eye of the subject
%cameraFs: frame rate of pupil tracking camera
%nFramesToRead: number of frames collected per trial
%WRITTEN BY: KYLE GHERES UPDATED: 12-08-2021 KWG

% This code requires a GPU for optimal performance
close all;
%% Check for existing mask of eye region
if ~exist('BW','var')
    BW=[];
    intensityThresh=[];
    thresh_ok='n';
else
    thresh_ok='y';
end
if ~exist('filNum','var')
    filNum=1;
end
%% Get Files to analyze
if ~ischar(pupilCamFileID)
pupilCamFileID = uigetfile('*_PupilCam.bin','MultiSelect','off'); %If the variable 'pupilCamFileID' is not a character string of a filename this will allow you to manually select a file
end
fid = fopen(pupilCamFileID); % This reads the binary file in to the work space
fseek(fid,0,'eof'); %find the end of the video frame
fileSize = ftell(fid); %calculate file size
fseek(fid,0,'bof'); %find the begining of video frames

%% Constants
medFilt_Params=[11 11]; % [x y] dimensions for 2d median filter of images
theangles=(1:1:180); %projection angles measured during radon transform of pupil
pupilHistEdges=[1:1:256]; %Camera data is unsigned 8bit integers. Ignore 0 values
RadonThresh=0.05; % Arbitrary threshold used to clean up radon transform above values ==1 below ==0
PupilThresh=0.35; % Arbitrary threshold used to clean up inverse radon transform above values ==1 below ==0
BlinkThresh=0.35;% Arbitrary threshold used to binarize data for blink detection above values ==1 below ==0
imageHeight=200; %How many pixels tall is the frame
imageWidth=200;%How many pixels wide is the frame
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame; 
if ~exist('nFramesToRead','var')
    nFramesToRead=9000; %Number of frames collected in imaging trial
end
if ~exist('cameraFs','var')
    cameraFs=30;
end
imageStack = zeros(imageHeight,imageWidth,nFramesToRead);

%% Read .bin File to imageStack
for a = 1:nFramesToRead
%     disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid,a*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,a) =imrotate(img,-90);
end
imageStack=uint8(imageStack);% convert double floating point data to unsignned 8bit integers
%% Empty Structures
Pupil_Area(1:size(imageStack,3))=NaN; %Area of pupil
Pupil_Major(1:size(imageStack,3))=NaN;% Length of major axis of pupil
Pupil_Minor(1:size(imageStack,3))=NaN;% Length of minor axis of pupil
Pupil_Centroid(1:size(imageStack,3),2)=NaN; %Center of pupil
PupilBoundary(1:size(imageStack,1),1:size(imageStack,2),1:size(imageStack,3))=NaN;
PupilPix=cell(1,size(imageStack,3));
%Overlay(1:size(imageStack,1),1:size(imageStack,2),(1:3),1:size(imageStack,3))=NaN; %RGB image stack with movie overlayed with Pupil location and size

%% Select ROI containing pupil
WorkingImg=imcomplement(imageStack(:,:,2)); %grab frame from image stack
if isempty(BW) %If no predefined ROI was fed to function
    fprintf('Draw roi around eye\n')
    figure(201);
    annotation('textbox',[0.4,0.9,0.1,0.1],'String','Draw ROI around eye','FitBoxToText','on','LineStyle','none','FontSize',16);
    [BW]=roipoly(WorkingImg);
end
%% Set Pupil intensity threshold    
Thresh_Set=4.5; % stardard deviations beyond mean intensity to binarize image for pupil tracking
WorkingImg=imcomplement(imageStack(:,:,2)); %grab frame from image stack
FiltImg=medfilt2(WorkingImg,medFilt_Params); %median filter image
ThreshImg=uint8(double(FiltImg).*BW); %Only look at pixel values in ROI


[phat,pci]=mle(reshape(ThreshImg(ThreshImg~=0),1,numel(ThreshImg(ThreshImg~=0))),'distribution','Normal'); %This models the distribution of pixel intensities as a gaussian...
%and is used to estimate and isolate the population of pixels that
%contains the pupil

figure(101);
pupilHist=histogram(ThreshImg((ThreshImg~=0)),'BinEdges',pupilHistEdges);
xlabel('Pixel intensities');
ylabel('Bin Counts');
title('Histogram of image pixel intensities')

normCounts=pupilHist.BinCounts./sum(pupilHist.BinCounts); %Normalizes bin count to total bin counts
theFit=pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); %Generate distribution from mle fit of data
normFit=theFit./sum(theFit); %Normalize fit so sum of gaussian ==1

if isempty(intensityThresh) %If no intensity threshold was fed to function calculate threshold to remove sclera
intensityThresh=phat(1)+( Thresh_Set*phat(2)); % set threshold as 4 sigma above population mean estimated from MLE
end

testImg=ThreshImg;
testImg(ThreshImg>=intensityThresh)=1;
testImg(ThreshImg<intensityThresh)=0;
testThresh=labeloverlay(imageStack(:,:,1),testImg);

figure(102);plot(pupilHist.BinEdges(2:end),normCounts,'k','LineWidth',1);
xlabel('Pixel intensities');
ylabel('Normalized bin counts');
title('Normalized histogram and MLE fit of histogram');
hold on;
plot(pupilHist.BinEdges,normFit,'r','LineWidth',2);
xline(intensityThresh,'--c','LineWidth',1);
legend({'Normalized Bin Counts','MLE fit of data','Pixel intensity threshold'},'Location','northwest');
xlim([0 256]);

figure(104);imshow(testThresh);
title('Pixels above threshold');

if strcmpi(thresh_ok,'n')
thresh_ok=input('Is pupil threshold value ok? (y/n)\n','s');
end
if strcmpi(thresh_ok,'n')
    intensityThresh=input('Manually set pupil intensity threshold\n');
    testImg(ThreshImg>=intensityThresh)=1;
    testImg(ThreshImg<intensityThresh)=0;
    testThresh=labeloverlay(imageStack(:,:,1),testImg);
    figure(103);imshow(testThresh);
    title('Pixels above threshold');
end

fprintf('Running\n')
procStart=tic;

%% Process camera frames
imageFrames=imageStack;
roiInt(1:size(imageFrames,3))=NaN;
for framenum=1:size(imageStack,3)
    FiltImg=medfilt2(imcomplement(imageFrames(:,:,framenum)),medFilt_Params);%medfilt2(imcomplement(imageStack(:,:,framenum)),[11 11]);
    ThreshImg=uint8(double(FiltImg).*BW); %Only look at pixel values in ROI
    roiInt_temp=sum(ThreshImg,1);
    roiInt(framenum)=sum(roiInt_temp,2);
    isoPupil=ThreshImg;
    isoPupil(isoPupil<intensityThresh)=0;
    isoPupil(isoPupil>=intensityThresh)=1;
    isoPupil=medfilt2(isoPupil,medFilt_Params);
    RadPupil=radon(isoPupil);
    minPupil=min(RadPupil,[],1);
    minMat=repmat(minPupil,size(RadPupil,1),1);
    MaxMat=repmat(max((RadPupil-minMat),[],1),size(RadPupil,1),1);
    NormPupil=(RadPupil-minMat)./MaxMat; %Normalize each projection angle to its min and max values. Each value should now be between [0 1]
    ThreshPupil=NormPupil;
    ThreshPupil(NormPupil>=RadonThresh)=1;
    ThreshPupil(NormPupil<RadonThresh)=0; %Binarize radon projection
    RadonPupil=iradon(double(ThreshPupil),theangles,'linear','Hamming',size(WorkingImg,2)); %transform back to image space
    [Pupil_Pix,Pupil_Boundary]=bwboundaries(RadonPupil>PupilThresh*max(RadonPupil(:)),8,'noholes'); %find area corresponding to pupil on binary image
    FillPupil=Pupil_Boundary;
    FillPupil=imfill(FillPupil,'holes'); %fill any subthreshold pixels inside the pupil boundary
    area_filled=regionprops(FillPupil,'FilledArea','Image','FilledImage','Centroid','MajorAxisLength','MinorAxisLength');
    if size(area_filled,1)>1
        for num=1:size(area_filled,1)
            theArea(num)=area_filled(num).FilledArea;
        end
        maxArea=max(theArea);
        areaLogical=theArea==maxArea;
        area_filled=area_filled(areaLogical);
    end
    if ~isempty(area_filled)
        Pupil_Area(framenum)=area_filled.FilledArea;
        Pupil_Major(framenum)=area_filled.MajorAxisLength;
        Pupil_Minor(framenum)=area_filled.MinorAxisLength;
        Pupil_Centroid(framenum,:)=area_filled.Centroid;
        PupilBoundary(:,:,framenum)=FillPupil;
        Hold=labeloverlay(imageStack(:,:,framenum),FillPupil,'Transparency',0.8);
        Overlay(:,:,:,framenum)=Hold;
    else
        Pupil_Area(framenum)=NaN;
        Pupil_Major(framenum)=NaN;
        Pupil_Minor(framenum)=NaN;
        Pupil_Centroid(framenum,:)=NaN;
        PupilBoundary(:,:,framenum)=FillPupil;
        Hold=labeloverlay(imageStack(:,:,framenum),FillPupil);
        Overlay(:,:,:,framenum)=repmat(Hold,1,1,3);
    end
end
PupilTracker.Pupil_Area=Pupil_Area;
PupilTracker.Pupil_Major=Pupil_Major;
PupilTracker.Pupil_Minor=Pupil_Minor;
PupilTracker.Pupil_Centroid=Pupil_Centroid;
PupilTracker.Pupil_Boundary=PupilBoundary;
PupilTracker.Pupil_Pix=PupilPix;
PupilTracker.Overlay=Overlay;
PupilTracker.Eye_ROI=BW;
PupilTracker.ROI_intensity=roiInt;
proceEnd=toc(procStart);
procMin=proceEnd/60;
minText=num2str(procMin);
procSec=round(str2double(minText(2:end))*60,0);
secText=num2str(procSec);
fprintf(['File processing time ' minText(1) ' min ' secText ' seconds\n'])
Blinks=find((abs(diff(PupilTracker.ROI_intensity))./PupilTracker.ROI_intensity(2:end))>=BlinkThresh)+1; %Returns a vector containing the indicies...

%% Determine when eyes are closed.
PlotPupilArea=PupilTracker.Pupil_Area;
PlotPupilArea(Blinks)=NaN;
BlinkTimes(1:size(PupilTracker.Pupil_Area,2))=NaN;
BlinkTimes(Blinks)=1;
BlinkTimes=BlinkTimes*(1.1*max(PupilTracker.Pupil_Area(:)));
PupilTime=(1:length(PupilTracker.Pupil_Area))/cameraFs;

%% Visualize pupil diameter and blink times
figure;plot(PupilTime,PlotPupilArea,'k','LineWidth',1);
hold on; scatter(PupilTime,BlinkTimes,50,'r','filled');
xlabel('Time (sec)');
ylabel('Pupil area (pixels)');
title('Pupil area changes');
legend('Pupil area','Eyes closed');
end