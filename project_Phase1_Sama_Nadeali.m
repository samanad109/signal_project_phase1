%Q1
%Part 1:
clc;
close all;
clear;
image=imread('D:\term\term 4\signal\project\phase 1\SS-Project-Phase1\q1.jpg');
k=1.3;

image1=im2double(image(:,:,1));
image1_fft=fft2(image1);
image1_fft_real=real(image1_fft(1:size(image1,1),1:size(image1,2)));
Hhp=1-1./(1+(10^6)./((abs(linspace(-1,1,size(image1,1))).^2)));
image1_fft_filtered=Hhp.*image1_fft;
image1_fft_sharpened=k.*image1_fft_filtered+1;
image1_sharpened=ifft2(image1_fft_sharpened);
image1_sharpened=real(image1_sharpened(1:size(image1,1),1:size(image1,2)));

image2=im2double(image(:,:,2));
image2_fft=fft2(image2);
image2_fft_real=real(image2_fft(1:size(image2,1),1:size(image2,2)));
Hhp=1-1./(1+(10^6)./((abs(linspace(-1,1,size(image2,1))).^2)));
image2_fft_filtered=Hhp.*image2_fft;
image2_fft_sharpened=k.*image2_fft_filtered+1;
image2_sharpened=ifft2(image2_fft_sharpened);
image2_sharpened=real(image2_sharpened(1:size(image2,1),1:size(image2,2)));

image3=im2double(image(:,:,3));
image3_fft=fft2(image3);
image3_fft_real=real(image3_fft(1:size(image3,1),1:size(image3,2)));
Hhp=1-1./(1+(10^6)./((abs(linspace(-1,1,size(image3,1))).^2)));
image3_fft_filtered=Hhp.*image3_fft;
image3_fft_sharpened=k.*image3_fft_filtered+1;
image3_sharpened=ifft2(image3_fft_sharpened);
image3_sharpened=real(image3_sharpened(1:size(image3,1),1:size(image3,2)));

image_sharpened=zeros(565,565,3);
image_sharpened(:,:,1)=image1_sharpened;
image_sharpened(:,:,2)=image2_sharpened;
image_sharpened(:,:,3)=image3_sharpened;

image_fft_real(:,:,1)=image1_fft_real;
image_fft_real(:,:,2)=image2_fft_real;
image_fft_real(:,:,3)=image3_fft_real;

imwrite(image_fft_real , 'q1-res1.jpg' , 'jpg');
imwrite(image_sharpened , 'q1-res2.jpg' , 'jpg');

%%
%part2
clc;
close all;
clear;
image=imread('D:\term\term 4\signal\project\phase 1\SS-Project-Phase1\q1.jpg');
k=0.3;

image1=image(:,:,1);
image1_fft=fftshift(fft2(double(image1)));
%calculate the high pass filter using the given formula:
%size(image1);  size(image1)=565 by 565
u=-282:1:282;
v=u;
[u,v]=meshgrid(u,v);
H=4.*pi.*pi.*(u.^2+v.^2).*image1_fft;
image1_final=abs(ifft2(fftshift(H)));

image2=image(:,:,2);
image2_fft=fftshift(fft2(double(image2)));
H=4*pi*pi*(u.^2+v.^2).*image2_fft;
image2_final=abs(ifft2(fftshift(H)));

image3=image(:,:,3);
image3_fft=fftshift(fft2(double(image3)));
H=4*pi*pi*(u.^2+v.^2).*image3_fft;
image3_final=abs(ifft2(fftshift(H)));

final_image=zeros(565,565,3);
final_image(:,:,1)=image1_final;
final_image(:,:,2)=image2_final;
final_image(:,:,3)=image3_final;

max_image=max(final_image,[],"all");
min_image=min(final_image,[],"all");
final_image=(final_image-min_image.*ones(565,565,3)).*(255./(max_image-min_image));
final_image=k.*final_image+1.*double(image);

max_image=max(final_image,[],"all");
min_image=min(final_image,[],"all");
final_image=(final_image-min_image.*ones(565,565,3)).*(255./(max_image-min_image));
final_image=uint8(final_image);

imwrite(final_image , 'q1-res3.jpg' , 'jpg');

%%
%Q2:
clc;
close all;
clear;
image=imread('D:\term\term 4\signal\project\phase 1\SS-Project-Phase1\pic.jpg');
image_adjusted=adapthisteq(image,"NumTiles",[2 2],"ClipLimit",0.012,"NBins",40000,"Range","full","Distribution","exponential","Alpha",0.1);
subplot(1,2,2);
imshow(image_adjusted);
subplot(1,2,1);
imshow(image);
imwrite(image_adjusted , 'q2_output.jpg' , 'jpg');

%%
%Q3:
clc;
close all;
clear;
einstein=imread('D:\term\term 4\signal\project\phase 1\SS-Project-Phase1\einstein.jpg');
marilyn=imread('D:\term\term 4\signal\project\phase 1\SS-Project-Phase1\marilyn.jpg');
size(einstein)
size(marilyn)
einstein=rgb2gray(einstein);
marilyn=marilyn(1:234,1:225);
einstein=imrotate(einstein, 10 , 'bilinear' , 'crop');
lp_filter=imgaussfilt(marilyn , 2);
hp_filter=1.3*einstein - imgaussfilt(einstein , 5);
result=1/2*lp_filter + 1/(2*1.3)*hp_filter;
subplot(2,3,1)
imshow(lp_filter);
subplot(2,3,2)
imshow(hp_filter);
subplot(2,3,4)
imshow(marilyn);
subplot(2,3,5)
imshow(einstein);
subplot(2,3,[3,6])
imshow(result);
imwrite(result , 'q3_output.jpg' , 'jpg');
