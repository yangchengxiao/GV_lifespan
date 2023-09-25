%% Calculate the GV, ROIV and GVtopo   YCX 2021.9.22
clear;clc;close all
sub_path='E:\GV_2021_9_16\lifespan\DATA\FunRawARWSCF';
mask_path='E:\Mask\BN_Atlas_246_3mm.nii';
% gv
result_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSCF\GV';
mkdir(result_path)
y_gv(sub_path,mask_path,result_path,3);

% Calculate the GV, ROIV and GVtopo   --- GSR
clear;clc;close all
sub_path='E:\GV_2021_9_16\lifespan\DATA\FunRawARWSCF';
mask_path='E:\Mask\BN_Atlas_246_3mm.nii';
% gv_gsr
result_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSCF\GV_GSR';
mkdir(result_path)
y_gv_gsr(sub_path,mask_path,result_path,3);


% Calculate the GS, ROIS and GStopo
clear;clc;close all
sub_path='E:\GV_2021_9_16\lifespan\DATA\FunRawARWSCF';
mask_path='E:\Mask\BN_Atlas_246_3mm.nii';
% gs
result_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSCF\GS';
mkdir(result_path)
y_gs(sub_path,mask_path,result_path,3);

%% power spectrum of GV
clear;clc;close all
GV_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSC\GV';
result_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSC\GV_PowerSpec';
mkdir(result_path)

TR=2;
fs=1/TR;
nfft= 2^nextpow2(225);
window=hamming(16); %海明窗
overlap=8; %数据重叠50%

y_PowerSpectrum(GV_path,result_path,window,overlap,nfft,fs)


%% Calculate the ReHo   YCX 2021.9.22
clear;clc;close all
sub_path='E:\GV_2021_9_16\lifespan\DATA\FunRawARWSCF';
mask_path='E:\Mask\BN_Atlas_246_3mm.nii';
% gv
result_path='E:\GV_2021_9_16\lifespan\result\FunRawARWSCF\GV';
mkdir(result_path)
y_gv(sub_path,mask_path,result_path,3);
