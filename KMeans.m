
img = imread('lena.png');
cd('..');
K = 3;
indices = linspace(1,size(img,2),K);
mu = reshape(img(uint8((size(img,1)+1)/2),round(indices),:),[K,3]);
mu = double(mu);
img_test = reshape(img,[],3);
img_test = double(img_test);
dist1 = zeros(size(img_test,1),K);
oldI = -1*ones(size(img_test,1),1);
display(mu)
for t=1:100
    s = 0;
    dist1 = zeros(size(img_test,1),K);
    for k=1:K
        dist1(:,k) = sqrt(sum((img_test - (ones(size(img_test,1),1)*(mu(k,:)))).^2,2));
    end
    [M,I] = min(dist1,[],2);
    s = sum(I~=oldI);
    oldI = I;
    for k=1:K
        mu(k,:) = sum(img_test((I == k),:),1)./size((I(I==k)),1);
    end
    if( s == 0)
        break;
    end
end

%Printing the image
cluster_val1 = reshape(I,[size(img,1), size(img,2)]);
R = cluster_val1;
G = cluster_val1;
B = cluster_val1;
img_out = cat(3,R,G,B);
for i=1:size(img_out,1)
    for j=1:size(img_out,2)
        k = img_out(i,j,1);
        img_out(i,j,1) = mu(k,1);
        img_out(i,j,2) = mu(k,2);
        img_out(i,j,3) = mu(k,3);
    end
end
img_int = uint8(img_out);
figure(1)

subplot(2,1,1)
imshow(img)
subplot(2,1,2)
imshow(img_int);