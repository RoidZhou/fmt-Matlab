close all;
clear;

%% Read Image
I1 = rgb2gray(imread('0001.jpg'));
I1 = imresize(I1, [480 480]);
I2 = imread('0020.jpg');

%% Generate Sample
roi = [151 250 251 400];
rot = 10 * pi / 180;
shift = [0 0];
scale = 1.24;
roic = [(roi(1)+roi(2))/2 (roi(3)+roi(4))/2];

% Crop Image Patch
I3 = I1(roi(1):roi(2), roi(3):roi(4), :);
I4 = zeros(size(I1,1), size(I1,2), size(I1,3));
I4(roi(1):roi(2), roi(3):roi(4), :) = I3;
I4 = uint8(I4);

% Rotation and Scale on center
A1 = [1 0 -roic(2);
      0 1 -roic(1);
      0 0 1];
A2 = [cos(rot) sin(rot) 0;
      -sin(rot) cos(rot) 0;
      0 0 1];
A3 = [scale 0 0;
      0 scale 0;
      0 0 1];
A4 = [1 0 roic(2);
      0 1 roic(1);
      0 0 1];

% Translate Matrix
A5 = [1 0 shift(1);
      0 1 shift(2);
      0 0 1];

% Affine Transformation
A = A1'*A2'*A3'*A4'*A5';
tform = affine2d(A);
I5 = imwarp(I4, tform);

% Generate Mask
Imask = zeros(size(I1,1), size(I1,2));
Imask(roi(1):roi(2), roi(3):roi(4)) = 1;
Imask = logical(Imask);
Imask = imwarp(Imask, tform);

% Crop Size
xx = logical(sum(Imask, 1));
yy = logical(sum(Imask, 2));
xxx = (xx * 1:length(xx));
yyy = (yy' * 1:length(yy));
xxx(xxx==0) = [];
yyy(yyy==0) = [];
roic2 = floor([mean(yyy) mean(xxx)]);
roi2 = [roic2(1)-size(I1,1)/2+1 roic2(1)+size(I1,1)/2 ...
        roic2(2)-size(I1,2)/2+1 roic2(2)+size(I1,2)/2];
roi2(1) = max(roi2(1),1);
roi2(2) = min(roi2(2),size(I5,1));
roi2(3) = max(roi2(3),1);
roi2(4) = min(roi2(4),size(I5,2));

% Crop Size
I6 = zeros(size(I1,1),size(I1,2),size(I1,3));
I6(1:(roi2(2)-roi2(1)+1),1:(roi2(4)-roi2(3)+1),:) = ...
    (I5(roi2(1):roi2(2), roi2(3):roi2(4), :));
I6 = uint8(I6);

% Show Result
figure, imshow(I4);
figure, imshow(I6);

%% Fourier-Mellin
%% Fourier Transform
%Ia = rgb2gray(I1);
%Ib = rgb2gray(I6);
Ia = I1;
Ib = I6;
sizeX = size(Ia, 1);
sizeY = size(Ia, 2);
sizeC = size(Ia, 3);

% Convert to FFT
FA = fftshift(fft2(Ia));
FB = fftshift(fft2(Ib));

% High Pass Filter
HF = hipass_filter(sizeX, sizeY, sizeC);
IA = HF.*abs(FA);
IB = HF.*abs(FB);


% Visualization
if sizeC == 3
figure, subplot(131), imagesc(IA(:,:,1));
subplot(132), imagesc(IA(:,:,2));
subplot(133), imagesc(IA(:,:,3));
figure, subplot(131), imagesc(IB(:,:,1));
subplot(132), imagesc(IB(:,:,2));
subplot(133), imagesc(IB(:,:,3));
else
figure, imagesc(IA);
figure, imagesc(IB);
end

%% Log Polar Transform
Ntheta = 360;
Rtheta = [0 2*pi];
Nrho = 500;
Rrho = [1 min(sizeX/2,sizeY/2)];
%Rrho = [1 sqrt((sizeX/2)^2+(sizeY/2)^2)];

% Caculate Polar Matrix
theta = linspace(Rtheta(1), Rtheta(2), Ntheta+1); theta(end)=[];
rho = logspace(log10(Rrho(1)), log10(Rrho(2)), Nrho);
xx = rho'*cos(theta) + sizeY/2;
yy = rho'*sin(theta) + sizeX/2;
mask= (xx>sizeY) | (xx<1) | (yy>sizeX) | (yy<1);

method = 'cubic'; % 'nearest','linear','cubic'
if sizeC == 3
    r = interp2(IA(:,:,1), xx, yy, method);
    g = interp2(IA(:,:,2), xx, yy, method);
    b = interp2(IA(:,:,3), xx, yy, method);
    r(mask) = 0;
    g(mask) = 0;
    b(mask) = 0;
    L1(:,:,3) = b;
    L1(:,:,2) = g;
    L1(:,:,1) = r;
    r = interp2(IB(:,:,1), xx, yy, method);
    g = interp2(IB(:,:,2), xx, yy, method);
    b = interp2(IB(:,:,3), xx, yy, method);
    r(mask) = 0;
    g(mask) = 0;
    b(mask) = 0;
    L2(:,:,3) = b;
    L2(:,:,2) = g;
    L2(:,:,1) = r;
else
    L1 = interp2(IA(:,:,1), xx, yy, method);
    L1(mask) = 0;
    L2 = interp2(IB(:,:,1), xx, yy, method);
    L2(mask) = 0;
end

%L1(1:400,:) = [];
%L2(1:400,:) = [];

% visualization
figure, imagesc(L1(:,:,1));
figure, imagesc(L2(:,:,1));

%% Phase Coralation
FL1 = fft2(L1);
FL2 = fft2(L2);

FCross = exp(1i*(angle(FL1)-angle(FL2)));
FPhase = real(ifft2(FCross));
R = sum(FPhase,3);
R(:,90:270) = 0;
R(100:400,:) = 0;
%R = FPhase(:,:,1);

[x, y] = find(R==max(max(R)));
x(2:end) = []; y(2:end)=[];
angle = theta(y) * 180 / pi;
if x>250
    x=501-x;
    scale = 1/rho(x);
else
    scale = rho(x);
end

% visualization
figure, mesh(R);

%% Recover Image
a = -angle * pi / 180;
s = 1/scale;
B1 = [cos(a) sin(a) 0;
      -sin(a) cos(a) 0;
      0 0 1];
B2 = [s 0 0;
      0 s 0;
      0 0 1];
tform = affine2d(B1'*B2');
II = imwarp(I6, tform);
figure, imshow(II);