%% Matlab code for the said task
% importing image
i = imread('cameraman.tif');
img = im2double(i);
%% Generating noisy versions of the image
% Uniform noise with mean=40 and standard deviation=20
m = 40;
std = 20;
noisyimg1 = img + std*randn(size(img)) + m;
% Gaussian Noise with mean=40 and standard deviation=20
var = std*2;
noisyimg2 = imnoise(img,'gaussian',m,var);
% Salt & Pepper noise with density=0.1
d = 0.1;
noisyimg3 = imnoise(img,'salt & pepper',d);
figure;
subplot(3,1,1)
imshow(noisyimg1)
title('Uniform Distribution Noise');
subplot(3,1,2)
imshow(noisyimg2)
title('Gaussian Noise');
subplot(3,1,3)
imshow(noisyimg3)
title('Salt&Pepper');
% upon adding noises, it is observed that the first two noises almost wipes
% out the image while the third noise keeps less effect on the image.
%% Filters
% Lowpass Filter with r0 being 50,30 and 20
% Saving the size of the input_image in pixels-
% M : no of rows (height of the image)
% N : no of columns (width of the image)
[M, N] = size(noisyimg3);
FT_noisyimg1 = fft2(double(noisyimg1));
FT_noisyimg2 = fft2(double(noisyimg2));
FT_noisyimg3 = fft2(double(noisyimg3));
% Assign Cut-off Frequency change 50,30,20
D0 = 50;
% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
% getting spatial coordinates
[V, U] = meshgrid(v, u);
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = double(D <= D0);
% Convolution between the Fourier Transformed
% image and the mask
Gif1 = H.*FT_noisyimg1;
Gif2 = H.*FT_noisyimg2;
Gif3 = H.*FT_noisyimg3;
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image1 = real(ifft2(double(Gif1)));
output_image2 = real(ifft2(double(Gif2)));
output_image3 = real(ifft2(double(Gif3)));
% Displaying Input Image and Output Image
figure;
subplot(3, 1, 1), imshow(output_image1, [ ]),
subplot(3, 1, 2), imshow(output_image2, [ ]),
subplot(3, 1, 3), imshow(output_image3, [ ]), title('Ideal Low Pass Filter');
% conclusion is that upon D0=50 the ILPF shows better results.
%% -----------------------------------------------------------------------
% Butterworth Filters
% Butterworth lowpass filter according to the formula
n = 2;
D0 = 50;
H = 1./(1 + (D./D0).^(2*n));
% Convolution between the Fourier Transformed 
% image and the mask
Gb1 = H.*FT_noisyimg1;
Gb2 = H.*FT_noisyimg2;
Gb3 = H.*FT_noisyimg3;

% Getting the resultant image by Inverse Fourier Transform 
% of the convoluted image using MATLAB library function  
% ifft2 (2D inverse fast fourier transform)   
output_imageb1 = real(ifft2(double(Gb1)));
output_imageb2 = real(ifft2(double(Gb2)));
output_imageb3 = real(ifft2(double(Gb3)));
% Displaying Input Image and Output Image 
figure;
subplot(3, 1, 1), imshow(output_imageb1, [ ]), 
subplot(3, 1, 2), imshow(output_imageb2, [ ]),
subplot(3, 1, 3), imshow(output_imageb3, [ ]), title('Butterworth Filtered');
% conclusion is that after applying butterworth filter with cutoff
% frequency 50 image filtering is better done while keeping the order at 2.
% ------------------------------------------------------------------------
%% gaussian filter 
std = 30;
% designing gaussian filter 
H = exp(-(D.^2)./((2.*std).^2));
% convolving gaussian filter in every case of noisy images
Gg1 = H.*FT_noisyimg1;
Gg2 = H.*FT_noisyimg2;
Gg3 = H.*FT_noisyimg3;
% getting output image back in time domain
output_imageg1 = real(ifft2(double(Gg1)));
output_imageg2 = real(ifft2(double(Gg2)));
output_imageg3 = real(ifft2(double(Gg3)));
% plotting filtered images in each case
figure;
subplot(3, 1, 1), imshow(output_imageg1, [ ]), 
subplot(3, 1, 2), imshow(output_imageg2, [ ]),
subplot(3, 1, 3), imshow(output_imageg3, [ ]), title('Gaussian Filter');
% by varying standard deviation, resutlant image gets blurr. so its better
% maintained at std=30.
%% Magnitude spectra of all filters
% Ideal low pass filter magnitude spectrum
% plotting frequency
f = 0:1:M-1;
figure;
subplot(3,1,1), plot(f,abs(Gif1)), title('Magnitude of ILPF applied on Uniform noisy image');
subplot(3,1,2), plot(f,abs(Gif2)), title('Magnitude of ILPF applied on Gaussian noisy image');
subplot(3,1,3), plot(f,abs(Gif3)), title('Magnitude of ILPF applied on Salt&Pepper noisy image');
% Butterworth filter magnitude spectrum
figure;
subplot(3,1,1), plot(f,abs(Gb1)), title('Magnitude of Butterworth applied on Uniform noisy image');
subplot(3,1,2), plot(f,abs(Gb2)), title('Magnitude of Butterworth applied on Gaussian noisy image');
subplot(3,1,3), plot(f,abs(Gb3)), title('Magnitude of Butterworth applied on Salt&Pepper noisy image');
% Gaussian Filter magnitude spectrum
figure;
subplot(3,1,1), plot(f,abs(Gg1)), title('Magnitude of Gaussian applied on Uniform noisy image');
subplot(3,1,2), plot(f,abs(Gg2)), title('Magnitude of Gaussian applied on Gaussian noisy image');
subplot(3,1,3), plot(f,abs(Gg3)), title('Magnitude of Gaussian applied on Salt&Pepper noisy image');
% Conclusion = It is observed that the first two noises dirupts the
% information in our test images on a large scale which is not able to
% retrieve it at a later stage after applying all filters. All these
% filters work better on std = 50 beyond that the image starts loosing
% sharpness. 
% ------------------------------------------------------------------------