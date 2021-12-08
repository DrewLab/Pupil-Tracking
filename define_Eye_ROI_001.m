function[Eye_ROI,Pixel_Threshold]=define_Eye_ROI_001(pupilCamFileID)
%% READ ME
%INPUTS
%pupilCamFileID: character string of .bin file of eye camera

%OUTPUTS
%Eye_ROI: Logical matrix for locating eye position
%Pixel_Threshold:Pixel intensity threshold for binarizing pixels within ROI

%This scripts allows user defined creation of ROIs used for pupil tracking
%in 'PupilTracker_QZ_GPU_001' and pixel intensity thresholds for image
%binarization
%WRITTEN: KYLE GHERES UPDATED: 08-09-2021 KWG
close all
%% Constants
imageHeight=200; %How many pixels tall is the frame
imageWidth=200;%How many pixels wide is the frame
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame; 
pupilHistEdges=[1:1:256];
BW=[];
Thresh_Set=4.5; % stardard deviations beyond mean intensity to binarize image for pupil tracking
medFilt_Params=[5 5]; % [x y] dimensions for 2d median filter of images

%% Get Files to analyze
if ~ischar(pupilCamFileID)
pupilCamFileID = uigetfile('*_PupilCam.bin','MultiSelect','off'); %If the variable 'pupilCamFileID' is not a character string of a filename this will allow you to manually select a file
end
fid = fopen(pupilCamFileID); % This reads the binary file in to the work space
fseek(fid,0,'eof'); %find the end of the video frame
fileSize = ftell(fid); %calculate file size
fseek(fid,0,'bof'); %find the begining of video frames

%% Read .bin File to imageStack
for a = 1:10 %No reason to need more frames %nFramesToRead
%     disp(['Creating image stack: (' num2str(a) '/' num2str(nFramesToRead) ')']); disp(' ')
    fseek(fid,a*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,a) = imrotate(img,-90);
end
imageStack=uint8(imageStack);% convert double floating point data to unsignned 8bit integers

%% Select ROI containing pupil
WorkingImg=imcomplement(imageStack(:,:,2)); %grab frame from image stack
if isempty(BW)
    fprintf('Draw roi around eye\n')
    figure(201);
    annotation('textbox',[0.4,0.9,0.1,0.1],'String','Draw ROI around eye','FitBoxToText','on','LineStyle','none','FontSize',16);
    [BW]=roipoly(WorkingImg);
end

Eye_ROI=BW;

%% Set Pupil intensity threshold    
WorkingImg=imcomplement(imageStack(:,:,2)); %grab frame from image stack
FiltImg=medfilt2(WorkingImg,medFilt_Params); %median filter image
ThreshImg=uint8(double(FiltImg).*BW); %Only look at pixel values in ROI

[phat,~]=mle(reshape(ThreshImg(ThreshImg~=0),1,numel(ThreshImg(ThreshImg~=0))),'distribution','Normal'); %This models the distribution of pixel intensities as a gaussian...
%and is used to estimate and remove the population of pixels that
%make up the sclera of the eye

figure(101);
pupilHist=histogram(ThreshImg((ThreshImg~=0)),'BinEdges',pupilHistEdges);
xlabel('Pixel intensities');
ylabel('Bin Counts');
title('Histogram of image pixel intensities')

normCounts=pupilHist.BinCounts./sum(pupilHist.BinCounts); %Normalizes bin count to total bin counts
theFit=pdf('normal',pupilHist.BinEdges,phat(1),phat(2)); %Generate distribution from mle fit of data
normFit=theFit./sum(theFit); %Normalize fit so sum of gaussian ==1


intensityThresh=phat(1)+( Thresh_Set*phat(2)); % set threshold as 4 sigma above population mean estimated from MLE
testImg=ThreshImg;
testImg(ThreshImg>=intensityThresh)=1;
testImg(ThreshImg<intensityThresh)=0;
testThresh=labeloverlay(imageStack(:,:,1),testImg);

figure(102);plot(pupilHist.BinEdges(2:end),normCounts,'k','LineWidth',1);
xlabel('Pixel intensities');
ylabel('Normalized bin counts');
title('Normalized histogram and MLE fit of histogram');
hold on;
plot(pupilHist.BinEdges,normFit,'--c','LineWidth',2);
xline(intensityThresh,'-r','LineWidth',2);
legend({'Normalized Bin Counts','MLE fit of data','Pixel intensity threshold'},'Location','northwest');
xlim([0 256]);

figure(103);imshow(testThresh);
title('Pixels above threshold');

thresh_ok=input('Is pupil threshold value ok? (y/n)\n','s');

if strcmpi(thresh_ok,'n')
    intensityThresh=input('Manually set pupil intensity threshold\n');
    testImg(ThreshImg>=intensityThresh)=1;
    testImg(ThreshImg<intensityThresh)=0;
    testThresh=labeloverlay(imageStack(:,:,1),testImg);
    figure(103);imshow(testThresh);
    title('Pixels above threshold');
end
Pixel_Threshold=intensityThresh;
end
