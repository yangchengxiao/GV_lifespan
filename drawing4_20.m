%% The average GVtopo of 4 runs in HCP
clear;
load('F:\Projects\GV\HCP\result\FunRawCRSF\ICC\ICC_GV.mat');
result_path='F:\Projects\GV\fig\BrainMap\HCP';
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
header_path='G:\DataBase\SALD\FunRawARWSCF\sub-031274\Filtered_4DVolume.nii';
for irun=1:4
   fname=fullfile(result_path,['HCP_GVtopoRun',num2str(irun)]);
    y_vec2nii(TOPO_mean(:,irun),mask_path,header_path,fname)
end


%% Age pattern of GV/GS
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

%% Brain map in HCP
save_path='F:\Projects\GV\Figure\HCP';
mkdir(save_path)
cd(save_path)

nii_dir=dir(data_path);
nii_dir(1:2)=[];
for i=1:length(nii_dir)
    nii_path=fullfile(data_path,nii_dir(i).name);
    
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
    
    
    
    
    
    
    
    
    
    
    
    
