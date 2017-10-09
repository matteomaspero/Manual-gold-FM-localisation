%% ----------------------------- LoA_MoACM --------------------------------
% Goal: Calculation of the accuracy in terms of inter-marker distance
% among corresponding FMs
% The procedure is described in the paper:
% This method is provided withouht any warraty, and is presented for
% reproducibility purposes. The code cannot be used for medical purposed
% without any further check. No responsability can be related to the
% developers from any misuse of the code.
% Reproduction and modification of the code is possible upon reference to
% the original code and the published paper and as long as the new code is
% released under the same creative commons license.
% The code has been developed by Matteo Maspero in 2017 at UMC Utrecht.
% For info feel free to contact:
% m.maspero@umcutrecht.nl/matteo.maspero.it@gmail.com

%% Housekeeping

clc; clear
close all;

Flag='Single'; % Chose 'Multiple' or 'Single' sequences

% load('./MatFile/Pos.mat')
if strcmp(Flag,'Multiple')
    load('./Data/Pos_T1.mat')
elseif strcmp(Flag,'Single')
    load('./Data/Pos.mat')
else
    error('The spellign of Flag is uncorrect!')
end
load('./Data/CoGCT.mat','CoGC*')

Detected_mrman1=3*ones(1,17);
Detected_mrman2=Detected_mrman1; Detected_mrman3=Detected_mrman1;
Detected_mrman4=Detected_mrman1; Detected_mrman5=Detected_mrman1;

%% Ord CT accordingly to the MR localisations!

CoGCT_ord{2}=CoGCT{2};
ssl=5;
Direct=3;
[~,Ord1]=sort(CoGCT{ssl}(:,Direct));
CoGCT_ord{ssl}=CoGCT{ssl}(Ord1,:);

%% Relative Distance

Ptno=1:17;
for ii=1:17
    RelDistCT_nocorr(:,ii)=sqrt(sum(abs(diff([CoGCT_ord{ii};CoGCT_ord{ii}(1,:)])).^2,2));
    RelDistMR1_nocorr(:,ii)=sqrt(sum(abs(diff([CoGMRman1_ord{ii};CoGMRman1_ord{ii}(1,:)])).^2,2));
    RelDistMR2_nocorr(:,ii)=sqrt(sum(abs(diff([CoGMRman2_ord{ii};CoGMRman2_ord{ii}(1,:)])).^2,2));
    RelDistMR3_nocorr(:,ii)=sqrt(sum(abs(diff([CoGMRman3_ord{ii};CoGMRman3_ord{ii}(1,:)])).^2,2));
    RelDistMR4_nocorr(:,ii)=sqrt(sum(abs(diff([CoGMRman4_ord{ii};CoGMRman4_ord{ii}(1,:)])).^2,2));
    RelDistMR5_nocorr(:,ii)=sqrt(sum(abs(diff([CoGMRman5_ord{ii};CoGMRman5_ord{ii}(1,:)])).^2,2));
end

%if on single sequence
if strcmp(Flag,'Single')
    RelDistCT_nocorr(:,5)=RelDistCT_nocorr([2 3 1],5);
end
RelMR1_CT=RelDistMR1_nocorr-RelDistCT_nocorr;
RelMR2_CT=RelDistMR2_nocorr-RelDistCT_nocorr;
RelMR3_CT=RelDistMR3_nocorr-RelDistCT_nocorr;
RelMR4_CT=RelDistMR4_nocorr-RelDistCT_nocorr;
RelMR5_CT=RelDistMR5_nocorr-RelDistCT_nocorr;
MaskExclud=logical(ones(3,17));

%Exclusion Criteria comes from LoA_LoACM.m Tand the figure 3 of the paper
if strcmp(Flag,'Multiple')
    disp('Exclusion imprecisely localted FM on multiple sequence')
    MaskExclud(1:2,4)=0;        % FM2 excluded
    MaskExclud([1,3],6)=0;      % FM1 
    MaskExclud(1:2,9)=0;        % FM2
    MaskExclud(1:3,14)=0;       % FM1
    MaskExclud(:,14)=0;       % Pt 14 excluded since out of the statistics (hip implant)

elseif strcmp(Flag,'Single')
    disp('Exclusion imprecisely localted FM on Single sequence')
    MaskExclud(:,4)=0;          % Pt4  since FM1 and FM2 need exclusion
    MaskExclud([1,3],6)=0;      % FM1
    MaskExclud([1,2],7)=0;      % FM2
    MaskExclud(:,9)=0;          % Pt9 since FM1 and FM2 need exclusion
    MaskExclud(:,14)=0; % Pt 14 excluded since out of the statistics
    MaskExclud(:,17)=0;         % Pt17 since FM1, 2 and 3 need exclusion
end

disp('                             Mean / Median +- STD [min;max]')

A=RelMR1_CT;
fprintf('Relative Distance CT vs MR 1: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(A(MaskExclud))),median(abs(A(MaskExclud))),std(abs(A(MaskExclud))),...
    min(abs(A(MaskExclud))),max(abs(A(MaskExclud))))
A=RelMR2_CT;
fprintf('Relative Distance CT vs MR 2: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(A(MaskExclud))),median(abs(A(MaskExclud))),std(abs(A(MaskExclud))),...
    min(abs(A(MaskExclud))),max(abs(A(MaskExclud))))
A=RelMR3_CT;
fprintf('Relative Distance CT vs MR 3: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(A(MaskExclud))),median(abs(A(MaskExclud))),std(abs(A(MaskExclud))),...
    min(abs(A(MaskExclud))),max(abs(A(MaskExclud))))
A=RelMR4_CT;
fprintf('Relative Distance CT vs MR 4: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(A(MaskExclud))),median(abs(A(MaskExclud))),std(abs(A(MaskExclud))),...
    min(abs(A(MaskExclud))),max(abs(A(MaskExclud))))
A=RelMR5_CT;
fprintf('Relative Distance CT vs MR 5: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(A(MaskExclud))),median(abs(A(MaskExclud))),std(abs(A(MaskExclud))),...
    min(abs(A(MaskExclud))),max(abs(A(MaskExclud))))
AS=[RelMR1_CT;RelMR2_CT;RelMR3_CT;RelMR4_CT;RelMR5_CT];

MaskExclud=[MaskExclud;MaskExclud;MaskExclud;MaskExclud;MaskExclud];
fprintf('Relative Distance CT vs MR All: %0.1f / %.1f +- %.1f1, [%.1f, %.1f] \n',...
    mean(abs(AS(MaskExclud))),median(abs(AS(MaskExclud))),std(abs(AS(MaskExclud))),...
    min(abs(AS(MaskExclud))),max(abs(AS(MaskExclud))))

