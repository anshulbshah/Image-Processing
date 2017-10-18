%% 
%Non Local Means Filtering
g_og = im2double(imread('krishna.png'));
f_og = im2double(imread('krishna_0_001.png'));
f = padarray(g_og,[15,15],'both');
g = padarray(f_og,[15,15],'both');
figure(1);
subplot(1,2,1);
imshow(g);
subplot(1,2,2);
imshow(f);
pixels = size(f,1)*size(f,2)*size(f,3);
%% 

i_max = size(f,1);
j_max = size(f,2);
W = 5;
Wsim = 3;
f_hat = zeros(size(g));
sigmasq = 0.1
for i = 15:i_max-15
    for j = 15:j_max-15
        patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
        filter = zeros(2*W+1);
        for k = i-W:i+W
            for l = j-W:j+W
                patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
                filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./sigmasq);
            end
        end
        filter = filter./sum(sum(filter));
        f_hat(i,j,:) = reshape(g(i-W:i+W,j-W:j+W,:),[],3)'*reshape(filter,[],1);
    end
end
MSE = reshape((f-f_hat),[],1)'*reshape((f-f_hat),[],1)/pixels;
PSNR = 10*log10(1/MSE);
%% 

%% 
figure(2)
imshow(f_hat);
figure(3)
subplot(1,3,1)
imshow(g)
subplot(1,3,2)
imshow(f)
subplot(1,3,3)
imshow(f_hat)
%% 
%Question 1
%a
%% 
i_max = size(f,1);
j_max = size(f,2);
W=5
Wsim=3
PSNRs = zeros(10,1)
pixels = size(f,1)*size(f,2)*size(f,3);
o = 1
for sigma=0.1:0.1:1
    f_hat = zeros(size(g));
    for i = 15:i_max-15
        for j = 15:j_max-15
            patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
            filter = zeros(2*W+1);
            for k = i-W:i+W
                for l = j-W:j+W
                    patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
                    filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./sigma^2);
                end
            end
            filter = filter./sum(sum(filter));
            f_hat(i,j,:) = reshape(g(i-W:i+W,j-W:j+W,:),[],3)'*reshape(filter,[],1);
        end
    end
    MSE = reshape((f-f_hat),[],1)'*reshape((f-f_hat),[],1)/pixels;
    PSNRs(o) = 10*log10(1/MSE);
    o=o+1
end

figure;
plot(PSNRs);
W=10
Wsim=3
PSNRs2 = zeros(10,1)
pixels = 276*276*3;
i_max = size(f,1);
j_max = size(f,2);
o = 1
filter = zeros(2*W+1);
f_hat = zeros(size(g));
for sigma=0.1:0.1:1
    for i = 15:i_max-15
        for j = 15:j_max-15
            patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
            for k = i-W:i+W
                for l = j-W:j+W
                    patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
                    filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./sigma^2);
                end
            end
            filter = filter./sum(sum(filter));
            f_hat(i,j,:) = reshape(g(i-W:i+W,j-W:j+W,:),[],3)'*reshape(filter,[],1);
        end
    end
    MSE = reshape((f-f_hat),[],1)'*reshape((f-f_hat),[],1)/pixels;
    PSNRs2(o) = 10*log10(1/MSE);
    o=o+1;
end        
  
MSE = reshape((f-f_hat),[],1)'*reshape((f-f_hat),[],1)/pixels;
PSNR_base = 10*log10(1/MSE);
%% 
%Comparision with Gaussian Filtering
%% 
size_o = 5;
PSNR_Gauss=zeros(10,1)
o=1
for sigma=0.1:0.1:1
     [x,y] = meshgrid(-size_o:1:size_o,-size_o:1:size_o);
     kernel = (1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2));
     kernel = kernel./(sum(sum(kernel)));
     f_hat = cat(3,convolution_operation(g(:,:,1),kernel),convolution_operation(g(:,:,2),kernel),convolution_operation(g(:,:,3),kernel));
     MSE = reshape((f-f_hat),[],1)'*reshape((f-f_hat),[],1)/pixels;
     PSNR_Gauss(o) = 10*log10(1/MSE);
     o=o+1;
end

%% 
figure;
sigma=0.1:0.1:1
plot(sigma,PSNRs);
hold on;
plot(sigma,PSNRs2);
hold on;
plot(sigma,PSNR_Gauss);
legend('Case1','Case2','Gauss');
%% 
%% 
W=5
Wsim=3
sigma_nlm=0.5;
sigma = 1;
i = 15+63;
j= 15+93;
filter = zeros(2*W+1);
patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
for k = i-W:i+W
    for l = j-W:j+W
        patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
        filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./sigma_nlm^2);
    end
end
filter = filter./sum(sum(filter));
[x,y] = meshgrid(-size_o:1:size_o,-size_o:1:size_o);
kernel = (1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2));
kernel = kernel./(sum(sum(kernel)));
figure;
subplot(2,1,1)
imshow(kernel,'InitialMagnification',1000);
subplot(2,1,2)
imshow(filter,'InitialMagnification',1000);

W=5
Wsim=3
sigma_nlm=0.5;
sigma = 1;
i = 15+77;
j= 15+118;
filter = zeros(2*W+1);
patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
for k = i-W:i+W
    for l = j-W:j+W
        patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
        filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./sigma_nlm^2);
    end
end
filter = filter./sum(sum(filter));
[x,y] = meshgrid(-size_o:1:size_o,-size_o:1:size_o);
kernel = (1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2));
kernel = kernel./(sum(sum(kernel)));
figure;
subplot(2,1,1)
imshow(kernel,'InitialMagnification',1000);
subplot(2,1,2)
imshow(filter,'InitialMagnification',1000);
%% 
%Question 4

W = 5;
Wsim = 3;
f_hat = zeros(size(g));
size_o=5;
for i = 15:i_max-15
    for j = 15:j_max-15
        patch = reshape(g(i-Wsim:i+Wsim,j-Wsim:j+Wsim,:),[],1)';
        filter = zeros(2*W+1);
        for k = i-W:i+W
            for l = j-W:j+W
                patch_temp = reshape(g(k-Wsim:k+Wsim,l-Wsim:l+Wsim,:),[],1)';
                filter(k-i+W+1,l-j+W+1) = exp(-(patch-patch_temp)*(patch - patch_temp)'./0.25);
            end
        end
        filter = filter./sum(sum(filter));
        f_hat(i,j,:) = reshape(g(i-W:i+W,j-W:j+W,:),[],3)'*reshape(filter,[],1);
    end
end
[x,y] = meshgrid(-size_o:1:size_o,-size_o:1:size_o);
kernel = (1/(2*pi*sigma^2))*exp(-(x.^2+y.^2)/(2*sigma^2));
kernel = kernel./(sum(sum(kernel)));
f_gauss = cat(3,convolution_operation(g(:,:,1),kernel),convolution_operation(g(:,:,2),kernel),convolution_operation(g(:,:,3),kernel));


figure;
subplot(1,3,1);
imshow(f_hat(15+77-W:15+77+W,15+118-W:15+118+W,:));
subplot(1,3,2);
imshow(g(15+77-W:15+77+W,15+118-W:15+118+W,:));
subplot(1,3,3);
imshow(f_gauss(15+77-W:15+77+W,15+118-W:15+118+W,:));

figure;
subplot(1,3,1);
imshow(f_hat(15+63-W:15+63+W,15+93-W:15+93+W,:));
subplot(1,3,2);
imshow(g(15+63-W:15+63+W,15+93-W:15+93+W,:));
subplot(1,3,3);
imshow(f_gauss(15+63-W:15+63+W,15+93-W:15+93+W,:));