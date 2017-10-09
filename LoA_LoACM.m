%% ----------------------------- LoA_MoACM --------------------------------
% Goal: Visualisation of the multi-observer (manual) gold fiducial marker
% localisation.
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

%% Initialisation Params

% 1234 to save the filkse as pdf
Gen.savePng=1234;

% where to save the fig
filedir='./Figure';

% Pt to be visualised
Pt=17;

%% Parameters for Visualization
addpath('../utils/')
Size_ftn=20; Name_ftn='Times New Roman';
set(0,'DefaultAxesFontSize',Size_ftn,'DefaultAxesFontName',Name_ftn,'DefaultTextFontname','Times New Roman');
NumCol=12;
ColorObs=distinguishable_colors(5,{'w','k'});
ColorDist=distinguishable_colors(NumCol);
Marke={'.','x','+','o'};
MarkeSiz={16,6,8,6};

%% Data opening for multiple sequences

load('./Data/Pos_T1.mat')
Ptno=1:numel(CoGMRman1_ord);

%% Calculation LoA, LoACM & agreement Position

PosCent=zeros(3,3,17);
CM_CT=zeros(17,3);
CM_MRman1=CM_CT;CM_MRman2=CM_CT;CM_MRman3=CM_CT;CM_MRman4=CM_CT;CM_MRman5=CM_CT;
for ss=1:numel(Ptno)
    CoGAll=zeros(3,3,5);
    CM_CT(ss,:)=mean(CoGCT{ss});
    CM_MRman1(ss,:)=mean(CoGMRman1_ord{ss}(1:size(CoGMRman1_ord{ss},1),:));
    CM_MRman2(ss,:)=mean(CoGMRman2_ord{ss}(1:size(CoGMRman2_ord{ss},1),:));
    CM_MRman3(ss,:)=mean(CoGMRman3_ord{ss}(1:size(CoGMRman3_ord{ss},1),:));
    CM_MRman4(ss,:)=mean(CoGMRman4_ord{ss}(1:size(CoGMRman4_ord{ss},1),:));
    CM_MRman5(ss,:)=mean(CoGMRman5_ord{ss}(1:size(CoGMRman5_ord{ss},1),:));
    CM_Agre(ss,:)=mean([CM_MRman1(ss,:);CM_MRman2(ss,:);CM_MRman3(ss,:);...
        CM_MRman4(ss,:);CM_MRman5(ss,:)]);
end

CMAll=cat(3,CM_MRman1,CM_MRman2,CM_MRman3,CM_MRman4,CM_MRman5);
Mat_CentTri=zeros(17,3,5);

for ssl=1:17
    CoGAll=zeros(3,3);
    CoGAll(1:size(CoGMRman1_ord{ss},1),:,1)=CoGMRman1_ord{ssl}(1:size(CoGMRman1_ord{ss},1),:);
    CoGAll(1:size(CoGMRman2_ord{ss},1),:,2)=CoGMRman2_ord{ssl}(1:size(CoGMRman2_ord{ss},1),:);
    CoGAll(1:size(CoGMRman3_ord{ss},1),:,3)=CoGMRman3_ord{ssl}(1:size(CoGMRman3_ord{ss},1),:);
    CoGAll(1:size(CoGMRman4_ord{ss},1),:,4)=CoGMRman4_ord{ssl}(1:size(CoGMRman4_ord{ss},1),:);
    CoGAll(1:size(CoGMRman5_ord{ss},1),:,5)=CoGMRman5_ord{ssl}(1:size(CoGMRman5_ord{ss},1),:);
    
    % Loop to check consistency of the dimensions
    jlls=0;
    Complete=[0];
    if size(CoGMRman1_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=1;
    end
    if size(CoGMRman2_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=2;
    end
    if size(CoGMRman3_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=3;
    end
    if size(CoGMRman4_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=4;
    end
    if size(CoGMRman5_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=5;
    end
    
    PosCent(:,:,ssl)=mean(CoGAll(:,:,Complete),3 );
    PosCentTri(:,:,ssl)=mean(sum(CoGAll(:,:,Complete),1)/3,3);
    for ll=1:3
        for ii=1:5
            Mat_CentTri(ssl,ll,ii)=abs(CMAll(ssl,ll,ii)-CM_Agre(ssl,ll));
            Dist1(ll,ii,ssl)=norm(CoGAll(ll,:,ii)-PosCent(ll,:,ssl));
            Dist_CentTri(ssl,ll,ii)=norm(CMAll(ssl,ll,ii)-CM_Agre(ssl,ll));
            Mat_CoGX(ssl,ll,ii)=abs(CoGAll(ll,1,ii)-PosCent(ll,1,ssl));
            Mat_CoGY(ssl,ll,ii)=abs(CoGAll(ll,2,ii)-PosCent(ll,2,ssl));
            Mat_CoGZ(ssl,ll,ii)=abs(CoGAll(ll,3,ii)-PosCent(ll,3,ssl));
        end
        rLoA(ll,ssl)=1.96*(std(Dist1(ll,:,ssl)));
        rLoA_CM(ll,ssl)=1.96*(std(Mat_CentTri(ssl,ll,:)));
        rLoA_X(ll,ssl)=1.96*(std(Mat_CoGX(ssl,ll,Complete)));
        rLoA_Y(ll,ssl)=1.96*(std(Mat_CoGY(ssl,ll,Complete)));
        rLoA_Z(ll,ssl)=1.96*(std(Mat_CoGZ(ssl,ll,Complete)));
    end
end

%% Visualization Fiducials in 3D with triangle MR/CT & CT high Inter-observer

figure('Name','Representation Inter-observer on multiple sequences')
colormap(ColorDist)

for ll=1:3
    view([-20 45]);
    axis equal;hold on
    if ll<=size(CoGMRman1_ord{Pt},1)
        plot3(CoGMRman1_ord{Pt}(ll,1),CoGMRman1_ord{Pt}(ll,2),CoGMRman1_ord{Pt}(ll,3),'Color',ColorObs(1,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 1')
    end
    if ll<=size(CoGMRman2_ord{Pt},1)
        plot3(CoGMRman2_ord{Pt}(ll,1),CoGMRman2_ord{Pt}(ll,2),CoGMRman2_ord{Pt}(ll,3),'Color',ColorObs(2,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 2')
    end
    if ll<=size(CoGMRman3_ord{Pt},1)
        plot3(CoGMRman3_ord{Pt}(ll,1),CoGMRman3_ord{Pt}(ll,2),CoGMRman3_ord{Pt}(ll,3),'Color',ColorObs(3,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 3')
    end
    if ll<=size(CoGMRman4_ord{Pt},1)
        plot3(CoGMRman4_ord{Pt}(ll,1),CoGMRman4_ord{Pt}(ll,2),CoGMRman4_ord{Pt}(ll,3),'Color',ColorObs(4,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 4')
    end
    if ll<=size(CoGMRman5_ord{Pt},1)
        plot3(CoGMRman5_ord{Pt}(ll,1),CoGMRman5_ord{Pt}(ll,2),CoGMRman5_ord{Pt}(ll,3),'Color',ColorObs(5,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 5')
    end
    plot3(PosCent(ll,1,Pt),PosCent(ll,2,Pt),PosCent(ll,3,Pt),'k',...
        'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Agreement')
    if ll==1
        s=legend('show'); box(s,'off')
        set(s,'Position',[ 0.3    0.78    0.0751    0.1273]);
    end
end
xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]'); % title(['Pt ',num2str(Pt)])
caxis([0,3])

% Center of Mass
annotation(gcf,'textbox',[0.11 0.765 0.137202377040826 0.14654282359762],...
    ... %s.Position + [0 s.Position(4) 0 0] ,...
    'String',{'*   FM 1','x   FM 2','+   FM 3','o  Center Mass'},...
    'FontSize',Size_ftn,'FontName',Name_ftn,'LineStyle','none',...
    'FitBoxToText','on');

plot3(CM_MRman1(Pt,1),CM_MRman1(Pt,2),CM_MRman1(Pt,3),'Color',ColorObs(1,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 1')
plot3(CM_MRman2(Pt,1),CM_MRman2(Pt,2),CM_MRman2(Pt,3),'Color',ColorObs(2,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 2')
plot3(CM_MRman3(Pt,1),CM_MRman3(Pt,2),CM_MRman3(Pt,3),'Color',ColorObs(3,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 3')
plot3(CM_MRman4(Pt,1),CM_MRman4(Pt,2),CM_MRman4(Pt,3),'Color',ColorObs(4,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 4')
plot3(CM_MRman5(Pt,1),CM_MRman5(Pt,2),CM_MRman5(Pt,3),'Color',ColorObs(5,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 5')
plot3(PosCentTri(1,1,Pt),PosCentTri(1,2,Pt),PosCentTri(1,3,Pt),'k',...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Agreement')
set(gcf,'Position',[1 1 1680 968])
grid on
if Gen.savePng==1234
    printpdf(gcf,[filedir,'/InterObs_Pt',num2str(Pt,'%i'),'_Multiple.pdf'])
end

%% variable and Remove
rLoA_X_T1=rLoA_X;rLoA_Y_T1=rLoA_Y; rLoA_Z_T1=rLoA_Z;
rLoA_CM_T1=rLoA_CM;
clear CoG*
clear rLoA_X rLoA_CM rLoA_Y rLoA_Z

%% Open Single Sequence

load('./Data/Pos.mat')

%% Calculation LoA, LoACM & agreement Position

PosCent=zeros(3,3,17);
CM_CT=zeros(17,3);
CM_MRman1=CM_CT;CM_MRman2=CM_CT;CM_MRman3=CM_CT;CM_MRman4=CM_CT;CM_MRman5=CM_CT;

for ss=1:numel(Ptno)
    CoGAll=zeros(3,3,5);
    CM_CT(ss,:)=mean(CoGCT{ss});
    CM_MRman1(ss,:)=mean(CoGMRman1_ord{ss}(1:size(CoGMRman1_ord{ss},1),:));
    CM_MRman2(ss,:)=mean(CoGMRman2_ord{ss}(1:size(CoGMRman2_ord{ss},1),:));
    CM_MRman3(ss,:)=mean(CoGMRman3_ord{ss}(1:size(CoGMRman3_ord{ss},1),:));
    CM_MRman4(ss,:)=mean(CoGMRman4_ord{ss}(1:size(CoGMRman4_ord{ss},1),:));
    CM_MRman5(ss,:)=mean(CoGMRman5_ord{ss}(1:size(CoGMRman5_ord{ss},1),:));
    CM_Agre(ss,:)=mean([CM_MRman1(ss,:);CM_MRman2(ss,:);CM_MRman3(ss,:);...
        CM_MRman4(ss,:);CM_MRman5(ss,:)]);
end

CMAll=cat(3,CM_MRman1,CM_MRman2,CM_MRman3,CM_MRman4,CM_MRman5);
Mat_CentTri=zeros(17,3,5);

for ssl=1:17
    CoGAll=zeros(3,3);
    CoGAll(1:size(CoGMRman1_ord{ss},1),:,1)=CoGMRman1_ord{ssl}(1:size(CoGMRman1_ord{ss},1),:);
    CoGAll(1:size(CoGMRman2_ord{ss},1),:,2)=CoGMRman2_ord{ssl}(1:size(CoGMRman2_ord{ss},1),:);
    CoGAll(1:size(CoGMRman3_ord{ss},1),:,3)=CoGMRman3_ord{ssl}(1:size(CoGMRman3_ord{ss},1),:);
    CoGAll(1:size(CoGMRman4_ord{ss},1),:,4)=CoGMRman4_ord{ssl}(1:size(CoGMRman4_ord{ss},1),:);
    CoGAll(1:size(CoGMRman5_ord{ss},1),:,5)=CoGMRman5_ord{ssl}(1:size(CoGMRman5_ord{ss},1),:);
    jlls=0;
    Complete=[0];
    if size(CoGMRman1_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=1;
    end
    if size(CoGMRman2_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=2;
    end
    if size(CoGMRman3_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=3;
    end
    if size(CoGMRman4_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=4;
    end
    if size(CoGMRman5_ord{ss},1)==3
        jlls=jlls+1;
        Complete(jlls)=5;
    end
    
    PosCent(:,:,ssl)=mean(CoGAll(:,:,Complete),3 );
    PosCentTri(:,:,ssl)=mean(sum(CoGAll(:,:,Complete),1)/3,3);
    for ll=1:3
        for ii=1:5
            Mat_CentTri(ssl,ll,ii)=abs(CMAll(ssl,ll,ii)-CM_Agre(ssl,ll));
            Dist1(ll,ii,ssl)=norm(CoGAll(ll,:,ii)-PosCent(ll,:,ssl));
            Dist_CentTri(ssl,ll,ii)=norm(CMAll(ssl,ll,ii)-CM_Agre(ssl,ll));
            Mat_CoGX(ssl,ll,ii)=abs(CoGAll(ll,1,ii)-PosCent(ll,1,ssl));
            Mat_CoGY(ssl,ll,ii)=abs(CoGAll(ll,2,ii)-PosCent(ll,2,ssl));
            Mat_CoGZ(ssl,ll,ii)=abs(CoGAll(ll,3,ii)-PosCent(ll,3,ssl));
        end
        rLoA(ll,ssl)=1.96*(std(Dist1(ll,:,ssl)));
        rLoA_CM(ll,ssl)=1.96*(std(Mat_CentTri(ssl,ll,:)));
        rLoA_X(ll,ssl)=1.96*(std(Mat_CoGX(ssl,ll,Complete)));
        rLoA_Y(ll,ssl)=1.96*(std(Mat_CoGY(ssl,ll,Complete)));
        rLoA_Z(ll,ssl)=1.96*(std(Mat_CoGZ(ssl,ll,Complete)));
    end
end

%% Visualization Fiducials in 3D with triangle MR/CT & CT high Inter-observer

figure('Name','Representation Inter-observer on single sequence')

for ll=1:3
    view([-20 45]);
    axis equal;hold on
    % after order
    if ll<=size(CoGMRman1_ord{Pt},1)
        plot3(CoGMRman1_ord{Pt}(ll,1),CoGMRman1_ord{Pt}(ll,2),CoGMRman1_ord{Pt}(ll,3),'Color',ColorObs(1,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 1')
    end
    if ll<=size(CoGMRman2_ord{Pt},1)
        plot3(CoGMRman2_ord{Pt}(ll,1),CoGMRman2_ord{Pt}(ll,2),CoGMRman2_ord{Pt}(ll,3),'Color',ColorObs(2,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 2')
    end
    if ll<=size(CoGMRman3_ord{Pt},1)
        plot3(CoGMRman3_ord{Pt}(ll,1),CoGMRman3_ord{Pt}(ll,2),CoGMRman3_ord{Pt}(ll,3),'Color',ColorObs(3,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 3')
    end
    if ll<=size(CoGMRman4_ord{Pt},1)
        plot3(CoGMRman4_ord{Pt}(ll,1),CoGMRman4_ord{Pt}(ll,2),CoGMRman4_ord{Pt}(ll,3),'Color',ColorObs(4,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 4')
    end
    if ll<=size(CoGMRman5_ord{Pt},1)
        plot3(CoGMRman5_ord{Pt}(ll,1),CoGMRman5_ord{Pt}(ll,2),CoGMRman5_ord{Pt}(ll,3),'Color',ColorObs(5,:),...
            'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Obs 5')
    end
    
    
    plot3(PosCent(ll,1,Pt),PosCent(ll,2,Pt),PosCent(ll,3,Pt),'k',...
        'Marker',Marke{ll},'MarkerSize',MarkeSiz{ll},'LineWidth',1.5,'DisplayName','Agreement')
    r=max(max([rLoA_X(ll,Pt) rLoA_Y(ll,Pt) rLoA_Z(ll,Pt)])) ;
    col=ceil(r/2.5*NumCol);
    if col>NumCol
        col=NumCol;
    end
    
    if ll==1
        s=legend('show'); box(s,'off')
        set(s,'Position',[ 0.3    0.78    0.0751    0.1273]);
    end
end
xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]'); % title(['Pt ',num2str(Pt)])
caxis([0,3])

% Center of Mass
annotation(gcf,'textbox',[0.11 0.765 0.137202377040826 0.14654282359762],...
    ... %s.Position + [0 s.Position(4) 0 0] ,...
    'String',{'*   FM 1','x   FM 2','+   FM 3','o  Center Mass'},...
    'FontSize',Size_ftn,'FontName',Name_ftn,'LineStyle','none',...
    'FitBoxToText','on');

plot3(CM_MRman1(Pt,1),CM_MRman1(Pt,2),CM_MRman1(Pt,3),'Color',ColorObs(1,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 1')
plot3(CM_MRman2(Pt,1),CM_MRman2(Pt,2),CM_MRman2(Pt,3),'Color',ColorObs(2,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 2')
plot3(CM_MRman3(Pt,1),CM_MRman3(Pt,2),CM_MRman3(Pt,3),'Color',ColorObs(3,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 3')
plot3(CM_MRman4(Pt,1),CM_MRman4(Pt,2),CM_MRman4(Pt,3),'Color',ColorObs(4,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 4')
plot3(CM_MRman5(Pt,1),CM_MRman5(Pt,2),CM_MRman5(Pt,3),'Color',ColorObs(5,:),...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Obs 5')
r=max(rLoA_CM(:,Pt));
col=ceil(r/2.5*NumCol);
if col>NumCol
    col=NumCol;
end
plot3(PosCentTri(1,1,Pt),PosCentTri(1,2,Pt),PosCentTri(1,3,Pt),'k',...
    'Marker','o','MarkerSize',6,'LineWidth',1.5,'DisplayName','Agreement')

set(gcf,'Position',[1 1 1680 968])

grid on
if Gen.savePng==1234
    printpdf(gcf,[filedir,'/InterObs_Pt',num2str(Pt,'%i'),'_Single.pdf'])
end

%% Single FM

th=2;
figure;set(gcf,'Position',[1           1        1680         960])
Col=distinguishable_colors(3);
subplot(321);
b=bar(rLoA_X');line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
title('Single sequence');
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
%legend({'FM 1' 'FM 2' 'FM 3','2 mm','1mm'},'Location','Best')%,'Position',[0.6818    0.8047    0.1486    0.1587]);
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
%legend('boxoff'); %title('Limit of Agreement');
ylabel('LoA_X [mm]','FontSize',Size_ftn,'FontName',Name_ftn)
xlim([0.5 17.5]);
set(gca,'Position',[0.05    0.67    0.455  0.28],'XTickLabel',{'','','','','','','',''})

subplot(322);
b=bar(rLoA_X_T1');line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
title('Multiple sequences');
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
legend({'FM 1' 'FM 2' 'FM 3','2 mm','1mm'},'Location','Best')%,'Position',[0.6818    0.8047    0.1486    0.1587]);
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
legend('boxoff'); %title('Limit of Agreement');
xlim([0.5 17.5]);
set(gca,'Position',[0.535    0.67    0.455  0.28],'XTickLabel',{'','','','','','','',''})

subplot(323); b=bar(rLoA_Y');ylabel('LoA_Y [mm]','FontSize',Size_ftn,'FontName',Name_ftn)
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
xlim([0.5 17.5]); line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
set(gca,'Position',[0.05    0.375    0.455  0.28],'XTickLabel',{'','','','','','','',''})

subplot(324); b=bar(rLoA_Y_T1');
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
xlim([0.5 17.5]); line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
set(gca,'Position',[0.535    0.375   0.455  0.28],'XTickLabel',{'','','','','','','',''})

subplot(325); b=bar(rLoA_Z'); ylabel('LoA_Z [mm]','FontSize',Size_ftn,'FontName',Name_ftn)
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
xlim([0.5 17.5]); line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
xlabel('Patient','FontSize',Size_ftn,'FontName',Name_ftn)
set(gca,'Position',[0.05    0.08    0.455  0.28])

subplot(326); b=bar(rLoA_Z_T1');
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.5 10.5])
xlim([0.5 17.5]); line([0.5 17.5],[th th],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th-1 th-1],'Color','k','LineStyle','--','LineWidth',1.5)
xlabel('Patient','FontSize',Size_ftn,'FontName',Name_ftn)
set(gca,'Position',[0.535 0.08  0.455  0.28])

if Gen.savePng==1234
    %saveas(gcf,[filedir,'/LoA_perFM_AllPt'],'png')
    %saveas(gcf,[filedir,'/LoA_perFM_AllPt.fig'])
    printpdf(gcf,[filedir,'/LoA_perFM_AllPt.pdf'])
    %  saveas(gcf,[filedir,'/LoA_perFM_AllPt.eps'],'epsc')
    % saveas(gcf,[filedir,'/LoA_perFM_AllPt.tif'])
end


%% CM

th=1;

figure;
set(gcf,'Position',[1           1        1640         500])
subplot(121)
b=bar(rLoA_CM');
title('Single sequence');
ylabel('LoA_{CM} [mm]'); xlabel('Patient')
xlim([0.5 17.5]);
line([0.5 17.5],[th+1 th+1],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th th],'Color','k','LineStyle','--','LineWidth',1.5)
%,'Position',[0.6818    0.8047    0.1486    0.1587]);
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.15 6])
%legend({'X' 'Y' 'Z','2 mm','1mm'},'Location','Best');
%legend('boxoff'); %title('Limit of Agreement');
ylabel('LoA_{CM} [mm]','FontSize',Size_ftn,'FontName',Name_ftn)
xlim([0.5 17.5]);
set(gca,'Position',[0.05 0.145 0.455  0.78])


subplot(122)
b=bar(rLoA_CM_T1');
title('Multiple sequences');
%ylabel('LoA_{CM} [mm]');
xlabel('Patient')
xlim([0.5 17.5]);
line([0.5 17.5],[th+1 th+1],'Color','k','LineStyle','-.','LineWidth',1.5)
line([0.5 17.5],[th th],'Color','k','LineStyle','--','LineWidth',1.5)
%,'Position',[0.6818    0.8047    0.1486    0.1587]);
for ii=1:3; b(ii).FaceColor=Col(ii,:); end
ylim([-0.15 6])
legend({'X' 'Y' 'Z','2 mm','1mm'},'Location','Best');
legend('boxoff'); %title('Limit of Agreement');
set(gca,'Position',[0.535 0.145 0.455  0.78])
%ylabel('X [mm]','FontSize',Size_ftn,'FontName',Name_ftn)
xlim([0.5 17.5]);

[a,b,c]=ttest(rLoA_CM_T1(:)',rLoA_CM(:)');

if Gen.savePng==1234
    %saveas(gcf,[filedir,'/LoA_CM_AllPt'],'png')
    %saveas(gcf,[filedir,'/LoA_CM_AllPt.fig'])
    printpdf(gcf,[filedir,'/LoA_CM_AllPt.pdf'])
    % saveas(gcf,[filedir,'/LoA_CM_AllPt.eps'],'epsc')
    % saveas(gcf,[filedir,'/LoA_CM_AllPt.tif'])
end

