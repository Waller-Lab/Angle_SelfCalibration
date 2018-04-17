% Fourier Ptychography Reconstruction 
% Main File
% Regina Eckert, 10/17/16
% Modified by Zack Phillips, 04/17/18

%Derived from code developed by Lei Tian, lei_tian@alum.mit.edu
%For ref: Lei Tian, et.al, Biomedical Optics Express 5, 2376-2389 (2014).

%% Input Parameters
addpath('Utilities');

%Input and output directories
base_dir = ['.' filesep 'data' filesep];

% This is the file header - _input or _output_selfCal will be appended
in_file = 'LED_cheekCell_comp_misaligned';

% Determine Calibration Type
% opts.position_calibration_type='no_calibration';
% opts.position_calibration_type='precalibrated';
% opts.position_calibration_type='online_nonrigid';
opts.position_calibration_type='online_rigid';
% opts.position_calibration_type='online_rigid_smooth';
% opts.position_calibration_type='simulated_annealing';
    
outDir=[base_dir 'results\'];
imgDir=base_dir;

% Determine if this method requires us to fit to a rigid transform or not
if contains(opts.position_calibration_type, 'online') && ~contains(opts.position_calibration_type, 'nonrigid')
    opts.fitRigidTransform = true;
else
    opts.fitRigidTransform = false;
end

newBk=false;
useUserNA=false;%true; %Use user-defined NA ***Takes precedence if both this and use_calibrationNA=true
use_calibrationNA=false; %Use calibrated NA
NA_user=0.26;

use_dpc_initialization=0;

% Options
nUsed=0; %0=all; -1=BF; >0, use that # sorted from lowest to highest NA
z0=0; %Defocus amount (um)
opts.tol = 1; %tol: maximum change of error allowed in two consecutive iterations
opts.maxIter = 10; %maxIter: maximum iterations 
opts.minIter = 2; %minIter: minimum iterations
opts.monotone = 1; %monotone (1, default): if monotone, error has to monotonically dropping when iters>minIter
opts.iters = 1; %Initialize iteration counter

% Display parameters
opts.display = 1; %display: display results (0: no (default) 1: yes)
opts.mode = 'real'; %mode: display in 'real' space or 'fourier' space

% Propagation type
opts.Prop_mode=0; %Propagation mode; if 0, Fresnel propagation; otherwise, angular spectrum

% Regularization & other processing controls
opts.OP_alpha = 10; %OP_alpha: regularization parameter for Object
opts.OP_beta =1; %OP_beta: regularization parameter for Pupil
opts.calbratetol = 1e-1; %calbratetol: parameter in controlling error tolerance in sa (simulated annealing)
opts.saveIterResult=false;

% Parameters for darkfield calibration
opts.scanRange=1;
opts.radialPenalty=0;
opts.gradIter=1;

% Parameters for removal of outliers (fitting to a rigid transform)
opts.transMethod='rigidScale';
opts.transAlpha=2;
opts.transScale=0.5;
opts.transTol=0.05;
opts.n_r=1;

%% FFT Operators
F = @(x) ifftshift(fft2(fftshift(x))); %Fourier Transform
Ft = @(x) ifftshift(ifft2(fftshift(x))); %Inverse Fourier Transform
upsamp = @(x) padarray(x,[(nObj(1)-NsampR(1))/2,(nObj(2)-NsampR(2))/2]);  

%% Load data
dataset=load([base_dir in_file '_input.mat']);

mag = dataset.metadata.objective.system_mag .* dataset.metadata.objective.mag;
NA_obj = dataset.metadata.objective.na;
dpix_c = dataset.metadata.camera.pixel_size_um;
lambda = dataset.metadata.illumination.wavelength_um;
zLED = dataset.metadata.illumination.z_distance_mm*10^3;

freqUV_design=dataset.metadata.source_list.na_design./lambda;
freqUV_cal=dataset.metadata.source_list.na_init./lambda;

NA_cal=dataset.metadata.self_cal.na_cal;
t_cal=dataset.metadata.self_cal.time_cal_s;
DFI=dataset.metadata.self_cal.DFI;

I=dataset.data;
NsampR=size(I);
numImg=NsampR(3);
NsampR=NsampR(1:2);

IbkThr=300;%2000; %Background threshold
if 1
    Bk= dataset.metadata.bk;
if newBk
    [bkC1, bkC2]=getRegions_bk(imgR(:,:,147)); %Displays sample image and asks user for input
    Bk(1,:) = mean(double(reshape(imgR(bkC1(1):bkC1(2),bkC1(3):bkC1(4),:),[],nImg))); 
    Bk(2,:) = mean(double(reshape(imgR(bkC2(1):bkC2(2),bkC2(3):bkC2(4),:),[],nImg)));
    
end
IBk= mean(Bk)'; %Take the mean across the two regions

if IBk(1)>IbkThr %Set the first one, if it's above the threshold
    IBk(1)=IbkThr;
end
for ii=2:length(IBk)
    if IBk(ii)>IbkThr
        IBk(ii)=IBk(ii-1); %If above the threshold, is brightfield; set to previous value
    end
end
else
    IBk = fpmData.metadata.bk;
end

I = I - repmat(permute(IBk,[3 2 1]),[NsampR 1]);
I(I<0)=0; %Where the region of interest < background, zero it out


%Make the output directory if it does not already exist
if ~exist(outDir, 'dir')
    mkdir(outDir) 
end


%% Load pre-calibration data
if contains(opts.position_calibration_type, 'precalibrate') || contains(opts.position_calibration_type, 'online') 
    pre_calibration_data = load([base_dir in_file '_output_selfCal.mat']);
    NA_cal=pre_calibration_data.metadata.self_cal.na_cal;
end

%% Process according to options

con=NsampR(1).*dpix_c./mag; %Conversion factor (pixels/(1/um))

%Figure out which illumination to use
if use_calibration
    freqUV=freqUV_cal;
else
    freqUV=freqUV_design;
end

illumNA=sqrt(sum(freqUV.^2,2)).*lambda; %abs(NA) of each illum
[illumNA,idx]=sort(illumNA); %Sort from smallest to largest
%Sort everything else as well
I=I(:,:,idx);
freqUV=freqUV(idx,:);
freqUV_cal=freqUV_cal(idx,:);
freqUV_design=freqUV_design(idx,:);
DFI=DFI(idx);

%Apply filter on which images to use
if nUsed==0
    %All
    idxUsed=1:numImg;
elseif nUsed<0
    %Brightfield
    idxUsed=~DFI;
else
    %Specified number
    idxUsed=1:nUsed;

end

%Select the images to be used
illumNA=illumNA(idxUsed);
I=I(:,:,idxUsed);
freqUV=freqUV(idxUsed,:);
DFI=DFI(idxUsed);
numImg=length(idxUsed); %Number of images

%Figure out which NA to use
if useUserNA 
    NA=NA_user; %User-defined
elseif use_calibrationNA
    NA=NA_cal; %Calibrated-defined
else
    NA=NA_obj; %System-defined
end

%% Set up system
[nObj, w_NA, x,y, u0, v0, k0]=defineSystem(lambda, NA, mag, dpix_c, NsampR, illumNA,opts.n_r);

sumI=sum(I,3)/size(I,3);


if use_dpc_initialization
    opts.O0 = initDPC( I, freqUV.*lambda, illumNA, NA_obj, lambda, dpix_c./mag, NsampR, nObj );
else
    opts.O0 = imresize(sqrt(sumI),nObj); %O0: initial guess for obj (in real space) %Initialize to normalized coherent sum
end

% Calibrate relative intensity of each LED
relInt=cos(asin(illumNA)) .^ 4.15; % Defined in Phillips et. al. 2015
opts.scale = ones(numImg,1); % LED brightness map (but all galvo presumed to be same value)
opts.P0 = w_NA; %P0: initial guess for P
opts.Ps = w_NA; %Pupil support
opts.H0 =defineProp(lambda, z0, k0, opts.Prop_mode); %Defocus kernel (or other known aberration)

opts.F = F; % F: operator of Fourier transform 
opts.Ft = Ft; % Ft: operator of inverse Fourier transform
opts.con=con; % Save the conversion factor

%spatial coordiates for object space
opts.x_obj=x;
opts.y_obj=y;

opts.dpix_mF=(dpix_c/mag)*NsampR(1)/nObj(1); %Pixel size of new image at sample plane
opts.freqUV_design=freqUV_design(idxUsed,:);

%% Reconstruct
[objR, Pupil, errorR, scaleR, freqUV_final, t_total] = fpmFunc(I, nObj, freqUV, opts);
%% Plot
param=outFile;
figH=figure(); imagesc(-angle(objR)); %caxis([-.6,0.6])
axis image; colormap gray; axis off; title(['Phase(Obj); param = ' param ]); colorbar; set(figH,'color','w')
% saveas(figH,[outDir outFile  '_phObjR'],'fig');
% export_fig([outDir outFile '_phObjR.png']);
figH=figure(); imagesc(abs(objR)); axis image; colormap gray; title(['abs(Obj); param = ' param ]); colorbar; set(figH,'color','w')
% saveas(figH,[outDir outFile '_absObjR'],'fig');
% export_fig([outDir outFile '_absObjR.png']);
% figH=figure(); imagesc(abs(F(objR))); axis image; title(['abs(F(Obj)); param = ' param ]); colorbar; caxis([0 1e5]); set(figH,'color','w')
% saveas(figH,[outDir outFile '_absObjF'],'fig');
% export_fig([outDir outFile '_absObjF.png']);
figH=figure(); imagesc(angle(Pupil)); axis image; title(['Phase(Pupil); param = ' param ]); colorbar; set(figH,'color','w'); colormap gray
figH=figure(); imagesc(abs(Pupil)); axis image; title(['Amplitude(Pupil); param = ' param ]); colorbar; set(figH,'color','w')
% % saveas(figH,[outDir outFile '_phPup'],'fig');
% % export_fig([outDir outFile '_phPup.png']);
% 
% [n,nbin]=hist(abs(objR(:).^2),100);
% [~,mI]=max(n);
% figH=figure(); imagesc((abs(objR.^2))./nbin(mI)); axis image; colormap gray; title(['Intensity; param = ' param ]); colorbar; set(figH,'color','w')
% caxis([0 2])
% saveas(figH,[outDir outFile '_int'],'fig');
% export_fig([outDir outFile '_int.png']);

figure(); plot(freqUV_design(:,1),freqUV_design(:,2),'*')
hold on; plot(freqUV_cal(:,1),freqUV_cal(:,2),'*')
hold on; plot(freqUV(:,1),freqUV(:,2),'*')
plot(freqUV_final(:,1),freqUV_final(:,2),'*')
axis image; legend('Designed','precalibratedibrated','Initialization','Final')%

%% Save

%Convert names to standard
obj=objR;
pupil=Pupil;
NA_used=NA;
%NA_obj
NA_self_cal=NA_cal;
cost=errorR;
dpix_recon=(dpix_c./mag)*NsampR(1)./nObj(1); %Pixel size in the reconstruction %pscrop
source_list_na=freqUV_final.*lambda;
source_list_na_design=freqUV_design.*lambda;
source_list_na_cal=freqUV_cal.*lambda;
source_list_na_init=freqUV.*lambda;
t_self_cal=t_cal;
%t_total
wavelength=lambda;
runVar.IbkThr=IbkThr;
runVar.OP_alpha=opts.OP_alpha;
runVar.OP_beta=opts.OP_beta;
runVar.maxIter=opts.maxIter;

%Check if outDir exists
if ~exist(outDir,'dir')
    mkdir(outDir)
end

%Save
save([outDir outFile '_results.mat'],'obj','pupil','NA_used','NA_obj','NA_self_cal','cost','dpix_recon',...
    'source_list_na','source_list_na_design','source_list_na_init','t_self_cal','t_total','wavelength','runVar')