%% Calculate the GV, ROIV and GVtopo   YCX 2021.9.22
clear;clc;close all
sub_path='G:\DataBase\SALD\FunRawARWSCF';
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
% gv
result_path='F:\Projects\GV\lifespan\result\FunRawARWSCF\GV';
mkdir(result_path)
y_gv(sub_path,mask_path,result_path,3);

% Calculate the GV, ROIV and GVtopo   --- GSR
clear;clc;close all
sub_path='G:\DataBase\SALD\FunRawARWSCF';
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
% gv_gsr
result_path='F:\Projects\GV\lifespan\result\FunRawARWSCF\GV_GSR';
mkdir(result_path)
y_gv_gsr(sub_path,mask_path,result_path,3);

% Calculate the GS, ROIS and GStopo
clear;clc;close all
sub_path='G:\DataBase\SALD\FunRawARWSCF';
mask_path='F:\Mask\BN_Atlas_246_3mm.nii';
% gs
result_path='F:\Projects\GV\lifespan\result\FunRawARWSCF\GS';
mkdir(result_path)
y_gs(sub_path,mask_path,result_path,3);

%% Calculate the GS, ROIS and GStopo  -- gvr
clear;clc;close all
gs_path='F:\Projects\GV\SALD\result\FunRawARWSCF\GS';
gv_path='F:\Projects\GV\SALD\result\FunRawARWSCF\GV';
result_path='F:\Projects\GV\SALD\result\FunRawARWSCF\GS_GVR';
mkdir(result_path);

sub_dir=dir(gs_path);
sub_dir(1:2)=[];
for isub=1:length(sub_dir)
    load(fullfile(gs_path,sub_dir(isub).name));   % load gs
    load(fullfile(gv_path,sub_dir(isub).name));   % load gv
    % regress
    x=[ones(length(gs),1),gv'];
    [~,~,gs]=regress(gs',x);
    gs=mapminmax(gs',min(gs),max(gs));    % normalization
    for iroi=1:size(rois,1)
        [~,~,rois(iroi,:)]=regress(rois(iroi,:)',x);
        rois(iroi,:)=mapminmax(rois(iroi,:)',min(rois(iroi,:)),max(rois(iroi,:)));    % normalization
        [r,p]=corrcoef(gs,rois(iroi,:));
        gs_topo.r(iroi)=r(2);
        gs_topo.p(iroi)=p(2);
    end
    gs_topo.z=a_fishertrans(gs_topo.r);
    fname=fullfile(result_path,sub_dir(isub).name);
    save(fname,'gs','rois','gs_topo')
end



