% 09.09.20
% starting point for BOLD amplitude analysis

clear all; close all; clc

addpath(genpath('Z:\COMMON')); % copy this folder to your machine and  point this path to it

%% load the data and show BOLD traces
%%%% preprocessed data in ROI space, with bandpassing
% suggest to start with bandpassed data given that that's what was used in
% the FC analysis
load('../data/precomputed/study1/fcAnalysisAfniProcAtlas_S1_bandpassed_data','roiEpochs','t1','t2','t3','subLabels','colors','TR','nSamples','nSubs','timeEpoch','fs','epochlen','epochstart',...
    'fcmatPre','fcmatDur','fcmatPost');

% data without bandpassing
% load('../data/precomputed/study1/fcAnalysisAfniProcAtlas_S1_data','roiEpochs','t1','t2','t3','subLabels','colors','TR','nSamples','nSubs','timeEpoch','fs','epochlen','epochstart',...
%     'fcmatPre','fcmatDur','fcmatPost');

% roiEpochs: echo time x subject x ROI x TR (time)
[nEchos,nSubjects,nRois,nTRs]=size(roiEpochs);

% each of the ROIs
% get the label of the first 5 ROIs
subLabels{1:5} % ROI is the "illuminated region" (estimated by MCML)
% ROIs 2-76 are in the left hemisphere, 77-151 are in the right hemisphere
% lh=left hemisphere, rh=right hemisphere
% G=gyrus
% S=sulcus

%% plot the time course of a sample echo/subject/roi:
e=2; s=1; r=1;
thisEpoch=squeeze(roiEpochs(e,s,r,:)); 
time=(0:nTRs-1)*TR/60;  % time in minutes
% NB: I cut out the first minute from the data, so the laser comes on at minute
% 9, not 10
figure; subplot(211); hold on
plot(time,thisEpoch); yl=ylim;
plot([9 9],[yl(1) yl(2)],'--r');
plot([19 19],[yl(1) yl(2)],'--k');
hlg=legend('BOLD','laser on','laser off');

%% look at the power spectrum of BOLD
nfft=2^nextpow2(nTRs);
freqs=(0:nfft-1)/nfft*(1/TR);
P=fft(roiEpochs,nfft,4); % fourth dimension is time
absP=abs(P);
tmp=permute(absP,[4 1 2 3]);
tmp=tmp(:,:);
muAbsP=nanmean(tmp,2); % grand average power spectrum of BOLD
figure;
subplot(221);
plot(freqs,muAbsP); xlim([0 0.5*1/TR]);

%% chow test analysis begins here
% pick a single echo/subject
thisRoiTs=squeeze(roiEpochs(1,1,:,:)); 
pvals=zeros(nRois,2);
tstats=zeros(nRois,2);
ALPHA=0.05;
CONSTANT=0; 
NLAGS=5; %p=NLAGS-1 
LASERONSETTR=find(time>9,1);
LASEROFFSETTR=find(time>19,1);
if CONSTANT % whether or not to add a constant term to AR model (?)
    coeffsPreStim=zeros(NLAGS,nRois,2); % intercept
    coeffsPrePost=zeros(NLAGS,nRois,2);
else
    coeffsPreStim=zeros(NLAGS-1,nRois,2); % no intercept
    coeffsPrePost=zeros(NLAGS-1,nRois,2);
end
for r=1:nRois
    r
    y=thisRoiTs(r,:);
    y=y(:); %needs to be a column vector
    X=tplitz(y,NLAGS); %  this creates a convolution matrix that is setting up the prediction of the current sample with the past
    if CONSTANT
        X=X(:,[end 2:end-1]); % include constant term (seems to cause numerical problems, so don't do it)
    else
        X=X(:,2:end-1); % don't include constant term 
    end
    
    %% pre versus stim
    [h1,p1,tstat1] = chowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASERONSETTR+1:LASEROFFSETTR,:); y2=y(LASERONSETTR+1:LASEROFFSETTR);
    coeffsPreStim(:,r,1)=X1\y1;
    coeffsPreStim(:,r,2)=X2\y2;
    
    %% pre versus post
    [h2,p2,tstat2] = chowtest(cat(1,X(1:LASERONSETTR,:),X(LASEROFFSETTR+1:end,:)),cat(1,y(1:LASERONSETTR),y(LASEROFFSETTR+1:end)),LASERONSETTR);
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASEROFFSETTR+1:end,:); y2=y(LASEROFFSETTR+1:end);
    coeffsPrePost(:,r,1)=X1\y1;
    coeffsPrePost(:,r,2)=X2\y2;
    
    pvals(r,1)=p1; pvals(r,2)=p2;
    tstats(r,1)=tstat1;  tstats(r,2)=tstat2; 
    
end
[pfdr,isSig]=fdr(pvals,ALPHA); % this is the correction for multiple comparisons (we are testing each roi individually)