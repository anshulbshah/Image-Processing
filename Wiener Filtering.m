%% 
sigma = 2
im = imread('lena.pgm');
im = im(2:end,2:end);
imshow(im);
blurred = imgaussfilt(im,sigma)
figure;
imshow(blurred)
J = im2double(imnoise(blurred,'gaussian',0,(1/255)^2));
figure;
imshow(J)
H = fft2(J);
F_2 = mat2gray(log(abs(H)));
imshow(F_2,[0,1]);
RMS_error = zeros(1991,1);
size1 = size(im,1);
sigma = 2;
size_o = 127;
[x,y] = meshgrid(-size_o:1:size_o,-size_o:1:size_o);
kernel = (1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2));
kernel = kernel./(sum(sum(kernel)));
G = fft2(kernel);

l=1;

%% 
figure;
imshow(G);
figure;
imshow(H);
figure;
imshow(J);
%% 
for k = 0.01:0.001:2

    Filtered = (conj(G)./(abs(G).^2 + k)).*H;
    F_2 = abs(ifft2(Filtered));
    %figure;
    F_2 = [F_2(size1/2:end,size1/2:end) F_2(size1/2:end,1:size1/2);
       F_2(1:size1/2,size1/2:end) F_2(1:size1/2,1:size1/2)];
    %imshow(F_2)
    RMS_error(l) = immse(F_2,im2double(im));
    l = l + 1;
end
[aa,bb] = min(RMS_error);
Filtered = (conj(G)./(abs(G).^2 + bb*0.001+0.01)).*H;
F_2 = ifft2(Filtered);
figure;
F_2 = [F_2(size1/2:end,size1/2:end) F_2(size1/2:end,1:size1/2);
   F_2(1:size1/2,size1/2:end) F_2(1:size1/2,1:size1/2)];
imshow(F_2)

figure;
imshow(J);