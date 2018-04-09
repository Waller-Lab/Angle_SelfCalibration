function [ freqDTh2 ] = calEstimate( FI, FIdiv, I, xI, yI, XYmid, radP, freqDTh, w_2NA, estTheta)
%Estimate the position of the circle centers based on thresholded image
%only

imSz=size(FI);
numImg=imSz(3);
imSz(3)=[];

centD=freqDTh(:,1);
theta=freqDTh(:,2);
freqXY=pol2Cart(freqDTh,XYmid);
XYopp=convOpp(freqXY, XYmid);


%Define pupil
pupil=sqrt((XYmid(1)-xI).^2 + (XYmid(2)-yI).^2)<=radP;
areaP=sum(pupil(:));

%Image processing
h=fspecial('average',3);
xLine=strel('line',3,0);
yLine=strel('line',3,90);
seD=strel('disk',2);

axB=round(imSz(1)/6); %How far around the x,y axes to look for lines
gradThr=3; %Threshold for gradient to find lines
yAx=false(imSz);
yAx(:,XYmid(1)-axB:XYmid(1)+axB)=1; %Map of y axis
xAx=false(imSz);
xAx(XYmid(2)-axB:XYmid(2)+axB,:)=1; %Map of x axis

%Initialization
xC2init=zeros(numImg,1);
yC2init=xC2init;
numC=5; %Number of images in imgC
imgSave1=zeros([imSz numImg]);
imgSave2=zeros([imSz numImg]);
thetaP=zeros(size(theta));

disp('1st Estimate:')
fprintf('Evaluating image: ')
strV=blanks(4);
fprintf('%s',strV) 
for ii=1:numImg
    
    %Print current image number
    str=num2str(ii);
    strV(1:length(str))=str;
    fprintf('\b\b\b\b%s',strV)
    
    %% Threshold log(F(I))
    %Threshold image
    lImg=log(abs(FI(:,:,ii)));
    mL=mean(lImg(:));
    sL=std(lImg(:));
    tL1=lImg>mL+sL;
    
    %**************************************************************************
    %Protect predicted area
    seed=(sqrt((freqXY(ii,1)-xI).^2+(freqXY(ii,2)-yI).^2)<=radP-0.165.*radP)|(sqrt((XYopp(ii,1)-xI).^2+(XYopp(ii,2)-yI).^2)<=radP-0.165.*radP);
    seed2=(sqrt((freqXY(ii,1)-xI).^2+(freqXY(ii,2)-yI).^2)<=radP)|(sqrt((XYopp(ii,1)-xI).^2+(XYopp(ii,2)-yI).^2)<=radP);
    %**************************************************************************
    
    %**************************************************************************
    %Will need tuning across different datasets/systems
    thMask=1+bwdist(seed).*0.000165.*radP;
    theL2=lImg>mL+thMask.*sL;
    %Remove the noise that's in the nonprotected area
    badPix=tL1&~theL2;
    %**************************************************************************

    %Smooth the image with an average filter
    %Empirically, shows same results as Weiner or median filter, but is
    %fastest
    mnImg=imfilter(lImg,h);
    mL=mean2(mnImg);
    sL=std2(mnImg);
    mnL1=mnImg>mL+sL; %Threshold
    mnL2=mnImg>mL+0.5.*sL;
    
    %Protect predicted area
    %**************************************************************************
    %Will need tuning across different datasets
    mtheL2=mnImg>mL+(thMask-0.5).*sL;
    badPixM=mnL2&~mtheL2;
    %**************************************************************************

    %% Threshold F(I)/|avg(F(I))|
    %Threshold the FI/avgFI image
    dImg=abs(FIdiv(:,:,ii));
    mD=mean(dImg(:));
    sD=std(dImg(:));
    tD1=dImg>mD+sD;

    %Threshold the mean FIdiv image
    mnDImg=imfilter(dImg,h);
    mD=mean2(mnDImg);
    sD=std2(mnDImg);
    mnD1=mnDImg>mD+sD;
    mnD2=mnDImg>mD+0.5.*sD;
    
    %Add the two together (often capture different areas)
    tB=logical(tL1+tD1);
    mnB=logical(mnL1+mnD1);

    %Combine into one dataset to be processed
    imgC=cat(3,tL1,mnL2,mnD2,tB,mnB);

    %% Remove noise
    %Initialize maps of "bad pixels" (lines at x & y axis)
    GxT=false(imSz); 
    GyT=GxT;
    for kk=1:size(imgC,3)
        %Find pixels with strong horiz, vert gradients in middle of image
        [Gx,Gy]=imgradientxy(imgC(:,:,kk));
        Gx=abs(Gx); Gy=abs(Gy); 

        %y-axis lines, Gx large
        Gxth=Gx>gradThr & yAx;
        Gxth2=imdilate(Gxth,xLine); %Dilate left-right to get middle
        GxT=GxT | Gxth2; %Combine into one map

        %x-axis lines, Gy large
        Gyth=Gy>gradThr & xAx;
        Gyth2=imdilate(Gyth,yLine);
        GyT=GyT | Gyth2;
    end

    basicBad=badPix | GxT | GyT | badPixM; %Combine maps of pixels to be discarded
    
    imgB=imgC & ~repmat(basicBad,[1 1 size(imgC,3)]); %Get rid of these pixels
    
    %% Clean up further
    %For image 3 (mean FIdiv, thresholded at 0.5*sigma) (mnD2)
    goodIm=imgB(:,:,3);
    goodIm=goodIm & ~bwareafilt(goodIm|seed2, [1 0.5.*areaP]) & w_2NA;
    goodIm=imfill(goodIm,'holes');
    
    goodIm2=imgB(:,:,5);
    goodIm2=goodIm2 & ~bwareafilt(goodIm2|seed2, [1 0.5.*areaP]) & w_2NA;
    goodIm2=imdilate(goodIm2,seD);
    goodIm2=imclose(goodIm2,seD);
    goodIm2=imfill(goodIm2,'holes');
    
    imgSave1(:,:,ii)=goodIm; %Save for further processing
    imgSave2(:,:,ii)=imgB(:,:,5); %Save for further processing
    
    %% Theta estimation
    if estTheta
        %Estimate theta
        sTemp=regionprops(goodIm2,'Orientation','MajorAxisLength','Eccentricity','Area');
        [~,aI]=max([sTemp.Area]); %Find max area
        angTemp=-sTemp(aI).Orientation; %Use only this blob
        if sTemp(aI).Area <= (areaP/2)
            %If less than half the area of one circle, likely center only,
            %which is perpendicular orientation to circles generally
            angTemp=angTemp+90;
        end

        %Make sure in same quadrant as original theta definition
        %***When we have a global shift, this becomes problematic for angles
        %that change dramatically (close to the center)
        if angTemp-theta(ii)>90
            thetaP(ii)=angTemp-180;
        elseif angTemp-theta(ii)<-90
            thetaP(ii)=angTemp+180;
        else
            thetaP(ii)=angTemp;
        end
    end
    

end
fprintf('\n')

%% Define values that go forward
%Estimating radius as radP
rad=radP;

%Decide whether to use estimated theta or the given theta when processing
if estTheta
    theta1=wrapTo180(thetaP(1:numImg));
else
    theta1=theta(1:numImg);
end



%% Estimate center distance

mIsave=zeros(numImg,1);
centD2init=mIsave;
guess=zeros(2,2,numImg);

%**************************************************************************
%Will need tuning across datasets and systems
dBnd=round([0.15.*imSz(1) 0.012.*imSz(1)]); %Initial bounds on the distance (lower, upper)
halfBnd=round(0.1281.*radP-3);
%**************************************************************************

%Get Radon transform size
temp=radon(abs(FI(:,:,ii)),0);
rSz=size(temp,1); %Length of Radon transform (varies with image size)
rCent=floor((rSz+1)/2); %Index of center
dCoor=(1:rSz)'; %Radon transform indices

disp('Center Distance Estimate:')
fprintf('Evaluating image: ')
strV=blanks(4);
fprintf('%s',strV) 
for ii=1:numImg
    
    %Print current image number
    str=num2str(ii);
    strV(1:length(str))=str;
    fprintf('\b\b\b\b%s',strV)
    
    use1=false;
    use2=false;
    
    %1st Estimate of Center via Radon Transform
    %Assume have the same theta, only estimate center distance
    %(Allow to vary in theta later)
    curCent=round(rCent-radP-centD(ii)); %Current center's predicted edge of circle (pixels) in Radon projection
    radT=radon(imgSave1(:,:,ii),180-theta1(ii)); %Project at 90 degrees to the predicted theta
    sradT=smooth(radT); %Smooth this
    dradT=diff(sradT); %Differentiate
    
    %Only accept points within -2 to +20 pixels of current center, and
    %where the function is increasing
    %***Values determined empirically for this system, and would likely have
    %to change for other systems
    %Try to find edge
    validLoc=dCoor>=curCent-dBnd(2) & dCoor<=curCent+dBnd(1) & [dradT; 0]>0;
    validLoc2=dCoor>=curCent-dBnd(2)+halfBnd & dCoor<=curCent+dBnd(1)+halfBnd;% & [dradT; 0]>0;
    
    %Find the side of the circle
    rt1=find(sradT>3 & validLoc,1); %Trying to find the edge; hardcoded as 3 because of noise
    rt2=find(sradT>(max(sradT)./2) & validLoc2,1); %Here trying to find point 4 pixels from the edge; amount found empirically, but makes sense because blurring with sigma=2;

    if isempty(rt1) || rt1+rad>rCent %Don't cross to the other side
        rt1=curCent;
        use2=true;
    end
    if isempty(rt2) || rt2+rad-halfBnd>rCent
        rt2=curCent;
        use1=true;
    end
    
    %Move in by a radius from the edge
    %Calculate x,y of center from these distance, theta pairs
    guess(1,:,ii)=[(rCent-(rt1+rad)).*cosd(theta1(ii))+XYmid(1),(rCent-(rt1+rad)).*sind(theta1(ii))+XYmid(2)];
    guess(2,:,ii)=[(rCent-(rt2+rad-halfBnd)).*cosd(theta1(ii))+XYmid(1),(rCent-(rt2+rad-halfBnd)).*sind(theta1(ii))+XYmid(2)];

    %Compare the RMSE error associated with each of these guesses
    imgBerr=-imageErr(I(:,:,ii),guess(:,1,ii),guess(:,2,ii),rad,xI,yI);
    [~,mI]=max(imgBerr); %Take the one with lower error (max of negative error)
    
    %But if we marked not to use something, change it
    %Will choose >3 + rad if there's a tie
    if use1 && mI==2
        mI=1;
    elseif use2 && mI==1
        mI=2;
    end
    
    %Save x,y coordinates (pixels)
    mIsave(ii)=mI;
    xC2init(ii)=guess(mI,1,ii);
    yC2init(ii)=guess(mI,2,ii);
    
    %Save center distance (pixels)
    guessCD=[rCent-(rt1+rad),rCent-(rt2+rad-halfBnd)];
    centD2init(ii)=guessCD(mI);
    
    centD2=centD2init;
    theta2=theta1;
    
    
end
freqDTh2=[centD2,theta2];

end

