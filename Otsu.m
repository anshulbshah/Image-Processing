%% 
clear all
close all

i1 = imread('palmleaf1.pgm');
i2 = imread('palmleaf2.pgm');

figure(1)
imshow(i1)

figure(2)
imshow(i2)

% figure(3)
% [level,EM] = graythresh(i2)
% BW = im2bw(i2,level);
% imshow(BW)
%% 
t_max = 0
max_val = -1000;
for i = 0:255
    [mu1,N1] = mu_1(i,i1);
    [mu2,N2] = mu_2(i,i1);
    [muT,NT] = mu_1(255,i1);
    sigb_sq = ((mu1-muT)^2)*N1/NT + ((mu2-muT)^2)*N2/NT;
    if sigb_sq > max_val
        max_val = sigb_sq;
        t_max = i;
    end
end
BW = im2bw(im2double(i1),t_max/255);
figure(3)
imshow(BW)


%% 
t_max = 0
max_val = -1000;
for i = 0:255
    [mu1,N1] = mu_1(i,i2);
    [mu2,N2] = mu_2(i,i2);
    [muT,NT] = mu_1(255,i2);
    sigb_sq = ((mu1-muT)^2)*N1/NT + ((mu2-muT)^2)*N2/NT;
    if sigb_sq > max_val
        max_val = sigb_sq;
        t_max = i;
    end
end
BW = im2bw(im2double(i2),t_max/255);
figure(4)
imshow(BW)