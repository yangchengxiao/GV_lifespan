%%  Check empty folder and delete them of dataset 'SouthwestSLIMData'
%%%  and find numbers of missing subject  2022.4.29    ycx
clc;close all;clear
data_path='E:\graduation_thesis\DATA\SouthwestSLIMData\Rest\Ses3\FunRawARW';
data_dir=dir(data_path);
data_dir(1:2)=[];

k=1;
for isub=1:length(data_dir)
    sub_dir=dir(fullfile(data_path,data_dir(isub).name));
    if length(sub_dir)==2
        empth_folder{k}=data_dir(isub).name;
        rmdir(fullfile(data_path,data_dir(isub).name));  % delete folder
        k=k+1;
    end
end
fname=fullfile('E:\graduation_thesis\DATA\SouthwestSLIMData\Rest\Ses3','DeleteFolder_empty.mat');
save(fname,'empth_folder');

%% Check three sessions integrity of each subject of dataset 'SouthwestSLIMData'
% and delete the imcomplete data   2022.5.25
clc;close all;clear
ses1_path='E:\GV_2021_9_16\SouthwestSLIMData\DATA\Rest\Ses1\FunRawARW';
ses2_path='E:\GV_2021_9_16\SouthwestSLIMData\DATA\Rest\Ses2\FunRawARW';
ses3_path='E:\GV_2021_9_16\SouthwestSLIMData\DATA\Rest\Ses3\FunRawARW';

ses1_dir=dir(ses1_path);
ses2_dir=dir(ses2_path);
ses3_dir=dir(ses3_path);

sub_num=1;
for i=3:length(ses1_dir)
    sub_name1=ses1_dir(i).name;
    for j=3:length(ses2_dir)
        sub_name2=ses2_dir(j).name;
        if strcmp(sub_name1,sub_name2)
            for k=3:length(ses3_dir)
                sub_name3=ses3_dir(k).name;
                if strcmp(sub_name2,sub_name3)
                   same_id{sub_num}=sub_name3;
                   sub_num=sub_num+1;
                end
            end
        end
    end
end

for i=3:length(ses1_dir)
    sub_name1=ses1_dir(i).name;
    if ~contains(sub_name1,same_id)
        rmdir(fullfile(ses1_path,sub_name1),'s'); 
    end
end

for j=3:length(ses2_dir)
    sub_name2=ses2_dir(j).name;
    if ~contains(sub_name2,same_id)
        rmdir(fullfile(ses2_path,sub_name2),'s'); 
    end
end

for k=3:length(ses3_dir)
    sub_name3=ses3_dir(k).name;
    if ~contains(sub_name3,same_id)
        rmdir(fullfile(ses3_path,sub_name3),'s'); 
    end
end

%% Move data to specified folder ---  NKI-RS   ycx 2022.4.7
clear;close all;clc
data_path='E:\DataBase\NKI-RS\RawData';
new_path='E:\GV_2021_9_16\NKI-RS\DATA\Rest1400_Raw1';
mkdir(new_path)

data_dir=dir(data_path);
data_dir(1:2)=[];

for i=1:length(data_dir)
    tic
    old_folder=fullfile(data_path,data_dir(i).name,'func');  % func TR=1400 ms
    if ~exist(old_folder,'dir')
        continue
    end
    
    func_dir=dir(old_folder);
    
    for ii=3:length(func_dir)
        if ~isempty(strfind(func_dir(ii).name,'rest_acq-1400_bold'))
            old_name=fullfile(old_folder,func_dir(ii).name);
            copyfile(old_name,new_path);
        end
    end
    a=strcat('copy: ',data_dir(i).name);
    disp(a);
    toc
end

%%% unzip and then
%%%
new_path2='E:\GV_2021_9_16\NKI-RS\DATA\Rest1400_Raw2';

data_dir=dir(new_path);
data_dir(1:2)=[];
wrong_label=1;
for isub=1:length(data_dir)
    tic
    a=strcat('copy1: ',num2str(isub));
    disp(a);
    old_name2=fullfile(data_dir(isub).folder,data_dir(isub).name);
    new_name2=fullfile(new_path2,data_dir(isub).name(1:13));
    try
        brain_dir=dir(old_name2);
        brain_name=fullfile(brain_dir(3).folder,brain_dir(3).name);
        [brain,~]=y_Read(brain_name);
        brain_size=size(brain);
        
        copyfile(old_name2,new_name2);
    catch
        sub_wrong{wrong_label}=data_dir(isub).name(1:13);
        wrong_label = wrong_label+1;
    end
    toc
end


%%  Check dataset 'NKI-RS' for completeness
%%%  and find numbers of missing subject  2022.4.29    ycx
clc;close all;clear
info_path='E:\graduation_thesis\DATA\NKI-RS\SubInfo_All.xlsx';
[~,~,sub_info]=xlsread(info_path);

data_path='I:\NKI-RS\RawDataBIDSLatest';
data_dir=dir(data_path);
data_dir(1:2)=[];

data_dir=struct2cell(data_dir);
data_name=data_dir(1,:)';

num_miss=1;
for isub=1:length(sub_info)
    if ~ismember(strcat('sub-',sub_info(isub)),data_name)
        sub_miss(num_miss,1)=sub_info(isub);
        num_miss=num_miss+1;
    end
end

%%  Check dataset 'NKI-RS' for integrity of .nii file and its .json fole(TR = 1400 / 645 ms)
%%%  and find numbers of missing subject  2022.4.29    ycx
clc;close all;clear    
TR='645';
info_path=fullfile('E:\GV_2021_9_16\NKI-RS\DATA',['SubInfo_REST',TR,'.xlsx']);
[~,~,sub_info]=xlsread(info_path);
sub_info(1,:)=[];
sub_names=sub_info(:,1);

data_path='I:\NKI-RS\RawDataBIDSLatest';
data_dir=dir(data_path);
data_dir(1:2)=[];

data_dir=struct2cell(data_dir);
data_name=data_dir(1,:)';

num_miss=1;
num_error1=1;   num_error2=1;
sub_error_nii=[];    sub_error_json=[];
for isub=1:length(sub_names)
    s=strcat('Processing: subnum - ',num2str(isub));
    disp(s);
    tic
    % Determine whether the subjects' data in the list has been downloaded, if not, record its number
    if ismember(strcat('sub-',sub_names(isub)),data_name)
        
        % Attempt to read the subject's .nii data, and record its number if it fails
        nii_name=strcat(strcat('sub-',sub_names(isub)),'_ses-BAS1_task-rest_acq-',TR,'_bold.nii.gz');
        sub_dir=char(fullfile(data_path,strcat('sub-',sub_names(isub)),'ses-BAS1','func',nii_name));
        try
            [brain,~]=y_Read(sub_dir);
            clear brain
            
        catch
            sub_error_nii(num_error1)=sub_info(isub);
            num_error1=num_error1+1;
        end
        
        % Attempt to read the subject's .json data, and record its number if it fails
        json_name=strcat(strcat('sub-',sub_names(isub)),'_ses-BAS1_task-rest_acq-',TR,'_bold.json');
        json_dir=char(fullfile(data_path,strcat('sub-',sub_names(isub)),'ses-BAS1','func',json_name));
        try
         data_json=loadjson(json_dir);   
         clear data_json   
        catch
            sub_error_json(num_error2)=sub_info(isub);
            num_error2=num_error2+1;
        end
    else
        sub_miss(num_miss,1)=sub_info(isub);
        num_miss=num_miss+1;
    end
    toc
end
fname=fullfile('E:\GV_2021_9_16\NKI-RS\DATA\miss_645','sub_miss.mat');
save(fname,'sub_error_json', 'sub_error_nii')


%%  copy file  --- NKI-RS  2022.7.4    ycx
clc;close all;clear    
TR='1400';
info_path=fullfile('E:\GV_2021_9_16\NKI-RS\DATA',['SubInfo_REST',TR,'.xlsx']);
[~,~,sub_info]=xlsread(info_path);
sub_info(1,:)=[];
sub_names=sub_info(:,1);

data_path='I:\NKI-RS\RawDataBIDSLatest';
data_dir=dir(data_path);
data_dir(1:2)=[];

data_dir=struct2cell(data_dir);
data_name=data_dir(1,:)';

out_path='E:\GV_2021_9_16\NKI-RS\DATA\REST_1400\FunRaw';
for isub=1:length(sub_names)
    s=strcat('Processing: subnum - ',num2str(isub));
    disp(s);
    tic
    % Determine whether the subjects' data in the list has been downloaded, if not, record its number
    if ismember(strcat('sub-',sub_names(isub)),data_name)
        
        out_dir=char(fullfile(out_path,strcat('sub-',sub_names(isub))));
        mkdir(out_dir);
        cd(out_dir);
        
        % Attempt to copy the subject's .nii data and .json data
        nii_name=strcat(strcat('sub-',sub_names(isub)),'_ses-BAS1_task-rest_acq-',TR,'_bold.nii.gz');
        nii_dir=char(fullfile(data_path,strcat('sub-',sub_names(isub)),'ses-BAS1','func',nii_name));
        copyfile(nii_dir);   % copy nii data
        
%         json_name=strcat(strcat('sub-',sub_names(isub)),'_ses-BAS1_task-rest_acq-',TR,'_bold.json');
%         json_dir=char(fullfile(data_path,strcat('sub-',sub_names(isub)),'ses-BAS1','func',json_name));
%         copyfile(json_dir);   % copy json data
    else
        sub_miss(num_miss,1)=sub_info(isub);
        num_miss=num_miss+1;
    end
    toc
end


%% Delete Files     --- NKI-RS
% The preprocessing steps are not all completed, 
% delete the completed subject's file
clc;close all;clear
data_path='I:\NKI_RS_prep\REST_1400\FunRaw';
data_dir=dir(data_path);
data_dir(1:2)=[];

for isub=1:length(data_dir)
    disp(num2str(isub))
    sub_dir=dir(fullfile(data_path,data_dir(isub).name));
    if length(sub_dir)==4
        if sub_dir(3).name(1) == "a"
        delete(fullfile(data_path,data_dir(isub).name,sub_dir(3).name));  % delete folder
        end
    end
end


%%  delete the small folder   2022.7.24
% Some subjects' fMRI data did not contain 404 time points, 
% which were deleted by the size of its NII file
clc;close all;clear
data_path='I:\NKI_RS_prep\REST_1400\FunRaw';
data_dir=dir(data_path);
data_dir(1:2)=[];
k=1;
for isub=1:length(data_dir)
    disp(num2str(isub))
    sub_dir=dir(fullfile(data_path,data_dir(isub).name));
    nii_size=sub_dir(3).bytes/(1000*1000);
    if nii_size < 1250
        wrong_folder{k,1}=data_dir(isub).name;
        wrong_folder{k,2}=nii_size;
        rmdir(fullfile(data_path,data_dir(isub).name),'s');  % delete folder
        k=k+1;
    end
end
fname=fullfile('I:\NKI_RS_prep\REST_1400','DeleteFolder_wrong.mat');
save(fname,'wrong_folder');


%%   data of 4 runs of HCP  is stored separately 2021.10.9
clc,clear;warning off all;close all
% read subject's data,
sub_path='E:\GV_2021_9_16\HCP\result\FunRawCRSF\GV_GSR\ALL';
save_path='E:\GV_2021_9_16\HCP\result\FunRawCRSF\GV_GSR\seperate';

save_path1=fullfile(save_path,'REST1_LR');
mkdir(save_path1)
save_path2=fullfile(save_path,'REST1_RL');
mkdir(save_path2)
save_path3=fullfile(save_path,'REST2_LR');
mkdir(save_path3)
save_path4=fullfile(save_path,'REST2_RL');
mkdir(save_path4)

cd(sub_path)
copyfile('*REST1_LR*',save_path1)
copyfile('*REST1_RL*',save_path2)
copyfile('*REST2_LR*',save_path3)
copyfile('*REST2_RL*',save_path4)


%% copy data of SALD to FunRaw
clc;clear
path_in='C:\Users\yangchengxiao\Downloads\Compressed\decompressed';
path_out='F:\DataBase\SALD\FunRaw';

sub_dir=dir(path_in);
sub_dir(1:2)=[];
for isub=124:length(sub_dir)
    isub
    data_dir=dir(fullfile(path_in,sub_dir(isub).name));
    
    copyname=fullfile(data_dir(4).folder,data_dir(4).name);
    out_path=fullfile(path_out,sub_dir(isub).name);
    mkdir(out_path)
    copyfile(copyname,out_path);
end




