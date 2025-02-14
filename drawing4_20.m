%% hcp 4 个run的平均GV，gvTOPO
clear;
load('F:\Projects\GV\HCP\result\FunRawCRSF\ICC\ICC_GV.mat');
result_path='F:\Projects\GV\fig\BrainMap\HCP';
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
header_path='G:\DataBase\SALD\FunRawARWSCF\sub-031274\Filtered_4DVolume.nii';
for irun=1:4
   fname=fullfile(result_path,['HCP_GVtopoRun',num2str(irun)]);
    y_vec2nii(TOPO_mean(:,irun),mask_path,header_path,fname)
end


%% GV / GS的年龄轨迹
clear;
load('F:\Projects\GV\lifespan\result\FunRawARWSCF\Age_effect\AgePattern_GV_GSR.mat');
figure(1);
subplot(1,2,1);
scatter(age,data_mean);
hold on
plot(fitresult_mean.power,'fit');
xlabel('Age');
ylabel('GV_m_e_a_n');
subplot(1,2,2)
scatter(age,data_std);
hold on
plot(fitresult_std.power,'fit');
xlabel('Age');
ylabel('GV_s_d');
set(gcf,'unit','centimeters','position',[20 10 22 8]);

[a,I]=min(feval(fitresult_mean.power,age));
age_min=age(I);


%% HCP 四个RUN的平均脑图-GVtopo/GStopo
clc;clear
result_path='F:\Projects\GV\HCP\result\FunRawCRSF\ICC';
cd(result_path)
data_path='F:\Projects\GV\HCP\result\FunRawCRSF\GS\seperate';
sub_path=dir(data_path);
sub_path(1:2)=[];
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
header_path='F:\DataBase\MDD-SPRBS\Data\SRPBS_MDD\KUT\Raw4DARWSCF\Band_0.0012_to_0.0728\sub-0670\CovRegressed_4DVolume_Filtered.nii';

for irun=1:4
    path_now=dir(fullfile(sub_path(irun).folder,sub_path(irun).name));
    path_now(1:2)=[];
    
    for isub=1:length(path_now)
        load(fullfile(path_now(isub).folder,path_now(isub).name));
        TOPO(isub,irun,:)=gs_topo.z;
    end
    TOPO_mean(irun,:)=mean(squeeze(TOPO(:,irun,:)),1);
    
    fname1=strcat('GS_',sub_path(irun).name);
    y_vec2nii(TOPO_mean(irun,:),mask_path,header_path,fname1)
end
[r_meantopo,p_meantopo]=corr(TOPO_mean');  % 平均TOPO的相关

% 所以被试四个RUN之间相关的平均
for isub=1:82
    [r,p]=corr(squeeze(TOPO(1,:,:))');
    R12(isub)=r(1,2);   % REST1_LR CORR REST1_RL
    R13(isub)=r(1,3);   % REST1_LR CORR REST2_LR
end
R12_mean=mean(R12);
R13_mean=mean(R13);

save CORR.mat TOPO_mean r_meantopo p_meantopo R12 R12_mean R13 R13_mean





%%   map Brain Map --- NKI-RS  --- each age group  2023.10.11
clc,clear;warning off all
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
header_path='F:\DataBase\MDD-SPRBS\Data\SRPBS_MDD\HRC\Raw4DARWSC\sub-0360\CovRegressed_4DVolume.nii';
data_path='F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_effect\AgePattern_GVtopo_GSR.mat';
load(data_path)
TOPO=TOPO_r';
age_key=[0,20,30,40,50,60,70,90];
output_path='F:\Projects\GV\fig\BrainMap\NKI-RS';
for ikey=1:7
    label=find(age>age_key(ikey) & age<=age_key(ikey+1));
    nii=mean(TOPO(label,:));
    % convert data to nii
    save_name=fullfile(output_path,strcat('GVCORR_GSR_group',num2str(ikey)));
    y_vec2nii(nii,mask_path,header_path,save_name)
end
%%% map nii of age effect --- NKI-RS  -- four kind age trend
Trend=["linear_pos";"linear_neg";"quard_pos";"quard_neg"];
Trend_label=[2,1,-2,-1];
for itrend=1:length(Trend)
    nii=zeros(1,246);
    nii(AgePattern_TOPO_r.choose==Trend_label(itrend))=AgePattern_TOPO_r.R2(AgePattern_TOPO_r.choose==Trend_label(itrend));
    % convert data to nii
    save_name=char(fullfile(output_path,strcat('GVCORR_GSR_agetrend_',Trend(itrend))));
    y_vec2nii(nii,mask_path,header_path,save_name)
end

%% brain map of TOPO's age pattern (the roi with max R)
clear;clc
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\AgePattern_GVTOPO.mat');
[~,key_gvtopo_pos]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-2));
gvtopo_pos=TOPO_r(key_gvtopo_pos,:)';
[~,key_gvtopo_neg]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-1));
gvtopo_neg=TOPO_r(key_gvtopo_neg,:)';
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\AgePattern_GVTOPO_GSR.mat');
% [~,key_gvtopo_gsr_pos]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-2));
gvtopo_gsr_pos=TOPO_r(key_gvtopo_pos,:)';
% [~,key_gvtopo_gsr_neg]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-1));
gvtopo_gsr_neg=TOPO_r(key_gvtopo_neg,:)';
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\AgePattern_GSTOPO.mat');
[~,key_gstopo_pos]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-2));
gstopo_pos=TOPO_r(key_gstopo_pos,:)';
[~,key_gstopo_neg]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-1));
gstopo_neg=TOPO_r(key_gstopo_neg,:)';
age=age';
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\AgePattern_GSTOPO_GVR.mat');
[~,key_gstopo_pos]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-2));
gstopo_pos=TOPO_r(245,:)';
[~,key_gstopo_neg]=max(AgePattern_TOPO_r.R2 .* (AgePattern_TOPO_r.choose==-1));
gstopo_neg=TOPO_r(197,:)';
age=age';

Out=y_Atlas246toAAL116([key_gvtopo_pos,key_gvtopo_neg,key_gstopo_pos,key_gstopo_neg]);

%%
save_path='F:\Projects\GV\Figure\HCP';
mkdir(save_path)
cd(save_path)

nii_dir=dir(data_path);
nii_dir(1:2)=[];
for i=1:length(nii_dir)
    nii_path=fullfile(data_path,nii_dir(i).name);
    %     data_dir=[data_path,filesep,['ICC_slow',num2str(i),'.nii']];
    
    [BrainNetViewerPath,fileN,extn]=fileparts(which('BrainNet.m'));
    SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
    CfgFile='F:\Projects\GV\Figure\HCP\option.mat';
    
    H_BrainNet=BrainNet_MapCfg(SurfFileName,nii_path,CfgFile);
    
    JpgFile=strcat('TOPO',nii_dir(i).name(1:end-4));
    %     JpgFile=strcat('ICC_slow_',num2str(i));
    %     eval(['print -r300 -djpeg -noui ''',JpgFile,''';']);
    eval(['print -r500 -dtiff -noui ''',JpgFile,''';']);
    
    close all
end
    
    
%% check the age pattern of GVtopo and GStopo
clear;clc
load('F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_effect\AgePattern_GVTOPO.mat');
AgePattern_gvtopo=AgePattern_TOPO_r.choose;
load('F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_effect\AgePattern_GSTOPO.mat');
AgePattern_gstopo=AgePattern_TOPO_r.choose;
label1=find(AgePattern_gvtopo==-2)';
label2=find(AgePattern_gstopo==-1)';

label=intersect(label1,label2);
Out=y_Atlas246toAAL116(label1);
    
    
    
    
    
    
    
    
    
    
    
    
    
