%%   ICC  --- HCP
%   Calculate the ICC (GV,GS) among 4 runs
%   2023.8.29
% ICC of GV/GS_mean, GV/GS_std
clc;clear;close all
data_path='F:\Projects\GV\HCP\result\FunRawCRSF\GV\seperate';
sub_path1=fullfile(data_path,'REST1_LR');
sub_path2=fullfile(data_path,'REST1_RL');
sub_path3=fullfile(data_path,'REST2_LR');
sub_path4=fullfile(data_path,'REST2_RL');
result_path='F:\Projects\GV\HCP\result\FunRawCRSF\ICC';
mkdir(result_path)
cd(result_path)
for irun=1:4
    data_dir=dir(eval(['sub_path',num2str(irun)]));
    data_dir(1:2)=[]; 
    sub_num=length(data_dir);
    for isub=1:sub_num
        sub_dir=fullfile(data_dir(isub).folder,data_dir(isub).name);
        load(sub_dir);
        DATA(:,irun,isub)=gv;  
        TOPO(:,irun,isub)=gv_topo.z;
    end
end
DATA_mean=squeeze(mean(DATA,1));
DATA_std=squeeze(std(DATA,1));
[icc_mean, LB, UB, F, df1, df2, p] = f_ICC(DATA_mean', '1-1');
[icc_std, LB, UB, F, df1, df2, p] = f_ICC(DATA_std', '1-1');

TOPO_mean=mean(TOPO,3);
[icc_topo2, LB, UB, F, df1, df2, p] = f_ICC(TOPO_mean, '1-1');
for iroi=1:246    % ICC of GVtopo, GStopo
    [icc_topo(iroi), LB, UB, F, df1, df2, p] = f_ICC(squeeze(TOPO(iroi,:,:))', '1-1');
end
% mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
% header_path='F:\DataBase\MDD-SPRBS\Data\SRPBS_MDD\HKH\Raw4DARWSCF\Band_0.0012_to_0.0728\sub-0425\CovRegressed_4DVolume_Filtered.nii';
% fname=fullfile(result_path,'ICC_GVtopo');
% y_vec2nii(icc_topo,mask_path,header_path,fname)

filename=fullfile(result_path,'ICC_GV');
save(filename,'icc_mean','icc_std','icc_topo','icc_topo2','TOPO_mean','DATA','DATA_mean','DATA_std')

%% the relationship between age and FD   ---NKI-RS  
clear;clc;close all
FD_path='G:\DataBase\SALD\RealignParameter';
FD_dir=dir(FD_path);
FD_dir(1:5)=[];
age_path='F:\Projects\GV\lifespan\sub_information.xlsx';  % load the information aboout age
info=xlsread(age_path);
%%% exclude sub whose mean head movement > 2 mm
[~,~,exclude_label]=xlsread('F:\Projects\GV\lifespan\data\exclude subjects_2mm.xlsx');
sub_num=1;
for isub=1:length(FD_dir)
    isub
    if ~sum(strcmp(FD_dir(isub).name,exclude_label))
        sub_path_fd=fullfile(FD_path,FD_dir(isub).name,['FD_Power_',FD_dir(isub).name,'.txt']);
        FD=load(sub_path_fd);
        FD_mean(sub_num)=mean(FD(2:end));
        
        label=find(string(num2str(info(:,1),'%05d')) == FD_dir(isub).name(end-4:end));
        age(sub_num)=info(label,3);   % sub's age
        gender(sub_num)=info(label,4);   % sub's gender
        sub_num=sub_num+1;
    else
        continue
    end
end
[r,p]=corrcoef(age,FD_mean);
Corr.r=r(2);
Corr.p=p(2);
result_path='F:\Projects\GV\lifespan\result\FunRawARWSCF\Age_effect';
mkdir(result_path);
cd(result_path);
save Corr_AgeAndFD.mat Corr FD_mean age gender


%% The age-effects of GV && GS  --- NKI-RS / lifespan
clear;clc
data_path='F:\Projects\GV\SALD\result\FunRawARWSCF\GS_GVR';
result_path='F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect';
mkdir(result_path)

% load the sub's age FD_mean
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\Corr_AgeAndFD.mat'); 

[~,~,exclude_label]=xlsread('F:\Projects\GV\SALD\data\exclude subjects_2mm.xlsx');
sub_num=1;

data_path_dir=dir(data_path);
data_path_dir(1:2)=[];
% data_path_dir(47)=[];   % NKI-RS loss gender information
for isub=1:length(data_path_dir)    % load gv of all subjects
    if ~sum(strcmp(data_path_dir(isub).name(1:end-4),exclude_label))
    load(fullfile(data_path,data_path_dir(isub).name))  
    data_mean(sub_num)=mean(gs) ;   % gv && gs
    sub_num=sub_num+1;
    else
        continue
    end
end
[GenderEffect_mean.h,GenderEffect_mean.p,~,GenderEffect_mean.stats] = ttest2(data_mean(gender==1),data_mean(gender==2));

%%% regress out FD_mean and gender
tbl = table(data_mean','VariableNames',"YY");   % mean
tbl.fd=FD_mean';
tbl.gender=categorical(gender');
mdl = fitlm(tbl,'YY ~  fd + gender');   % fit
data_mean_r=mdl.Residuals{:,4};   % residual
data_mean_r=mapminmax(data_mean_r',min(data_mean),max(data_mean));    % normalization
clear tbl
%%%% Fit 
FitType=["poly2","poly1"];   % poly2 is quadratic; poly1 is linear
for ifit=1:2
    ft=fittype(FitType(ifit));
    for idata=1:2
        switch idata
            case 1    % data_mean
                fit_data=data_mean';
            case 2   % data_mean_r
                fit_data=data_mean_r';
        end
        [power,gof]=fit(age', fit_data, ft);
        R2(ifit,idata)=gof.adjrsquare;
        F_value(ifit,idata)=y_AdjRsquare2F (gof.adjrsquare,2,length(age));
        P_value(ifit,idata)=1-fcdf(F_value(ifit,idata),2,length(age));
    end  
end
cd(result_path)
save AgePattern_GS_GVR.mat R2 F_value P_value data_mean data_mean_r  age


%% The age-effects of GVtopo && GStopo  --- NKI-RS / lifespan
clear;clc
data_path='F:\Projects\GV\SALD\result\FunRawARWSCF\GS_GVR';
result_path='F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect';
mkdir(result_path)
cd(result_path)
save_fname=fullfile(result_path,'AgePattern_GSTOPO_GVR.mat');
% load the sub's age FD_mean
load('F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect\Corr_AgeAndFD.mat'); 

[~,~,exclude_label]=xlsread('F:\Projects\GV\SALD\data\exclude subjects_2mm.xlsx');
sub_num=1;

data_path_dir=dir(data_path);
data_path_dir(1:2)=[];
% data_path_dir(47)=[];   % NKI-RS loss gender information
for isub=1:length(data_path_dir)
    if ~sum(strcmp(data_path_dir(isub).name(1:end-4),exclude_label))
    load(fullfile(data_path,data_path_dir(isub).name));
    TOPO(:,sub_num)=gs_topo.z;    % gvtopo  && gstopo
    sub_num=sub_num+1;
    else
        continue
    end
end
%%% regress out FD_mean and gender
for iroi=1:246
    tbl = table(TOPO(iroi,:)','VariableNames',"YY");   % mean
    tbl.fd=FD_mean';
    tbl.gender=categorical(gender');
    mdl = fitlm(tbl,'YY ~  fd + gender');   % fit
    TOPO_r(iroi,:)=mdl.Residuals{:,4};   % residual
    TOPO_r(iroi,:)=mapminmax(TOPO_r(iroi,:),min(TOPO(iroi,:)),max(TOPO(iroi,:)))';    % normalization
    clear tbl
end
%%% fit
FitType=["poly2","poly1"];   % poly2 is quadratic; poly1 is linear
for ifit=1:2
    ft=fittype(FitType(ifit));
    for idata=1:2
        switch idata
            case 1    % TOPO
                fit_data=TOPO;
            case 2    % TOPO_R
                fit_data=TOPO_r;
        end
        for iroi=1:246
            disp(num2str(iroi))
            [power,gof]=fit(age', fit_data(iroi,:)', ft);
            Pr(ifit,idata,iroi)=power.p1;
            R2(ifit,idata,iroi)=gof.adjrsquare;
            F_value(ifit,idata,iroi)=y_AdjRsquare2F (gof.adjrsquare,2,length(age));
            P_value(ifit,idata,iroi)=1-fcdf(F_value(ifit,idata,iroi),2,length(age));
        end
        % fdr
        P_value_fdr(ifit,idata,:)=mafdr(squeeze(P_value(ifit,idata,:)),'BHFDR','true');
    end
end
%
data_label=["TOPO","TOPO_r"];
for idata=1:2
    for iroi=1:246
        p= P_value_fdr(:,idata,iroi);
        if sum(p <= 0.05)  == 2
            [~,fit_label]=max(R2(:,idata,iroi));
        elseif  sum(p <= 0.05)  == 1
            [~,fit_label]=min(p);
        else
            fitresult_topo.R2(iroi)=0;
            fitresult_topo.p(iroi)=1;
            fitresult_topo.choose(iroi)=0;
            continue
        end
        fitresult_topo.R2(iroi)=R2(fit_label,idata,iroi);
        fitresult_topo.p(iroi)=P_value(fit_label,idata,iroi);
        if fit_label==1   %% quadratic
            if Pr(fit_label,idata,iroi) > 0
                fitresult_topo.choose(iroi)=-2;   % positive quadratic  -- 开口向上
            else
                fitresult_topo.choose(iroi)=-1;   % nagetive quadratic  - 开口向下
            end
        elseif fit_label==2   %% linear
            if Pr(fit_label,idata,iroi) > 0
                fitresult_topo.choose(iroi)=-2;   % positive linear
            else
                fitresult_topo.choose(iroi)=-1;   % nagetive linear
            end
        end
    end
    expr=['AgePattern_',char(data_label(idata)),'= ','fitresult_topo;';];
    eval(expr);   % 
end
save(save_fname,'AgePattern_TOPO','AgePattern_TOPO_r','TOPO','TOPO_r','age')

%% age-prediction   --- internal raliability
clear;clc
result_path='F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_Predict';
mkdir(result_path)
cd(result_path)
%%% load data
parent_path='F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_effect';
Folder=["AgePattern_GS.mat","AgePattern_GV.mat",...
    "AgePattern_GSTOPO.mat","AgePattern_GVTOPO.mat"];
Data1=[];   % save gs & gv 
for ifolder=1:2
    load(fullfile(parent_path,char(Folder(ifolder))));
    Data1=[Data1;data_mean_r];
end
Data2=[];   % save gs_topo & gv_topo 
for ifolder=3:4
    load(fullfile(parent_path,char(Folder(ifolder))));
    label_linear{ifolder-2}=find(AgePattern_TOPO_r.choose==1 | AgePattern_TOPO_r.choose==2);
    label_quard{ifolder-2}=find(AgePattern_TOPO_r.choose==-1 | AgePattern_TOPO_r.choose==-2);
    label_mix{ifolder-2}=find(AgePattern_TOPO_r.choose~=0 );
    Data2=[Data2;TOPO_r(label_mix{ifolder-2},:)];
end
% Create input matrix X and output vector y
y = age;
num_follds=10;
rounds=100;
numPermutations=5000;
Mat=[Data1;Data2];
K1=1;   % GSmean
K2=2;   % GVeman
K3=3:length(label_mix{1})+2;   % GSTOPO
K4=length(label_mix{1})+3 : length(label_mix{1})+length(label_mix{2})+2;   % GVTOPO

label_key{1}=[K1]; % GSmean
label_key{2}=[K2]; % GVmean
label_key{3}=[K1 K2]; % GVmean & GSmean
label_key{4}=[K1 K3]; % GVmean & GVTOPO
label_key{5}=[K2 K4]; % GSmean & GSTOPO
label_key{6}=[K2 K3 K4]; % GVmean & GSTOPO & GVTOPO 
label_key{7}=[K1 K2 K3 K4]; % GSmean & GVmean & GSTOPO & GVTOPO 
%%
for ilabel=1:7
    ilabel
    X=Mat(label_key{ilabel},:);
    [R(ilabel),MAE(ilabel)]=y_svm(X',y,'polynomial',num_follds,rounds);
%     [permCorrelations(ilabel,:),pValue(ilabel)]=y_PermutationTest(X',y,'polynomial',numPermutations,R(ilabel));
end
save AgePred_internal.mat R MAE  Mat label_linear label_quard label_mix label_key age

%%   age prediction     --- external reliability
clc;clear
%%% load train data --- NKI-RS
load('F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_Predict\AgePred_internal.mat');
Ytrain=age; clear age
%%% load test data --- SALD
parent_path='F:\Projects\GV\SALD\result\FunRawARWSCF\Age_effect';
Folder=["AgePattern_GS.mat","AgePattern_GV.mat",...
    "AgePattern_GSTOPO.mat","AgePattern_GVTOPO.mat"];
Data1=[];   % save gs & gv & gv_gsr
for ifolder=1:2
    load(fullfile(parent_path,char(Folder(ifolder))));
    Data1=[Data1;data_mean_r];
end
Data2=[];   % save gs_topo & gv_topo & gv_topo_gsr
for ifolder=3:4
    load(fullfile(parent_path,char(Folder(ifolder))));
    Data2=[Data2;TOPO_r(label_mix{ifolder-2},:)];
end
Xtest=[Data1;Data2];
Ytest = age;

for ilabel=1:7
    ilabel
    Xtrain=Mat(label_key{ilabel},:);
    svmModel = fitrsvm(Xtrain', Ytrain', 'KernelFunction', 'polynomial', 'PolynomialOrder', 2);
    ypred = predict(svmModel, Xtest(label_key{ilabel},:)');
    % Compute performance metrics
    r(ilabel) = corr(ypred, Ytest');
    mae(ilabel) = mean(abs(ypred - Ytest'));
    Ytest=Ytest';
end
result_path='F:\Projects\GV\NKI-RS\result\FunRawARWSCF\Age_Predict';
cd(result_path)
save AgePred_external.mat r mae



