close all
clear all
%Loading The images
i1 = rgb2gray(imresize(imread('ex2.jpg'),1/12));
i2 = rgb2gray(imresize(imread('ex1.jpg'),1/12));
%i3 = imread('img3.pgm');
figure(1)
imshow(i1);
figure(2)
imshow(i2);
%figure(3)
%imshow(i3);
%% 

%Determining the Homographies

[corresp1,corresp2] = sift_corresp('img1.pgm','img2.pgm');

figure(2)
subplot(1,2,1)
imshow(i1)
hold on;
plot(corresp1(:,2)',corresp1(:,1)','.');
subplot(1,2,2)
imshow(i2)
hold on;
plot(corresp2(:,2)',corresp2(:,1)','.');

no_of_iterations = 300;
dist_thresh = sqrt(10);
%Swapping Coords, so first column is x and second is y
% corresp1 = [corresp1(:,2)-size(i1,2)/2 corresp1(:,1)-size(i1,1)/2]; 
% corresp2 = [corresp2(:,2)-size(i2,2)/2 corresp2(:,1)-size(i2,2)/2];
corresp1 = [corresp1(:,2) corresp1(:,1)]; 
corresp2 = [corresp2(:,2) corresp2(:,1)];
for i=1:no_of_iterations
    ind = randsample(size(corresp1,1),4);
    A = [corresp2(ind(1),1) corresp2(ind(1),2) 1 0 0 0 -corresp1(ind(1),1)*corresp2(ind(1),1) -corresp1(ind(1),1)*corresp2(ind(1),2) -corresp1(ind(1),1);
         0 0 0 corresp2(ind(1),1) corresp2(ind(1),2) 1 -corresp1(ind(1),2)*corresp2(ind(1),1) -corresp2(ind(1),2)*corresp1(ind(1),2) -corresp1(ind(1),2);
         corresp2(ind(2),1) corresp2(ind(2),2) 1 0 0 0 -corresp1(ind(2),1)*corresp2(ind(2),1) -corresp1(ind(2),1)*corresp2(ind(2),2) -corresp1(ind(2),1);
         0 0 0 corresp2(ind(2),1) corresp2(ind(2),2) 1 -corresp1(ind(2),2)*corresp2(ind(2),1) -corresp2(ind(2),2)*corresp1(ind(2),2) -corresp1(ind(2),2);
         corresp2(ind(3),1) corresp2(ind(3),2) 1 0 0 0 -corresp1(ind(3),1)*corresp2(ind(3),1) -corresp1(ind(3),1)*corresp2(ind(3),2) -corresp1(ind(3),1);
         0 0 0 corresp2(ind(3),1) corresp2(ind(3),2) 1 -corresp1(ind(3),2)*corresp2(ind(3),1) -corresp2(ind(3),2)*corresp1(ind(3),2) -corresp1(ind(3),2);
         corresp2(ind(4),1) corresp2(ind(4),2) 1 0 0 0 -corresp1(ind(4),1)*corresp2(ind(4),1) -corresp1(ind(4),1)*corresp2(ind(4),2) -corresp1(ind(4),1);
         0 0 0 corresp2(ind(4),1) corresp2(ind(4),2) 1 -corresp1(ind(4),2)*corresp2(ind(4),1) -corresp2(ind(4),2)*corresp1(ind(4),2) -corresp1(ind(4),2)];
     H = reshape(null(A),[3,3])';
     [corresp1it,corresp2it] = deal(corresp1,corresp2);
     corresp1it(ind,:) = [];
     corresp2it(ind,:) = [];
     [x_rem,y_rem] = deal(corresp2it(:,1),corresp2it(:,2));
     [x_ac,y_ac] = deal(corresp1it(:,1)',corresp1it(:,2)');
     calculated = H*[x_rem';y_rem';ones(1,size(x_rem',2))];
     [x_calc,y_calc] = deal(calculated(1,:)./calculated(3,:),calculated(2,:)./calculated(3,:));
     dist = sqrt((x_calc-x_ac).^2+(y_calc-y_ac).^2);
     no_of_inliers = sum(dist<dist_thresh);
     if(no_of_inliers > 0.8*size(x_rem,1))
        break;
     end
end
figure(3)
imshow(i1);
hold on
plot(x_ac,y_ac,'.');
hold on
plot(x_calc,y_calc,'*');

H21 = H;
%% 
%% 
[corresp1,corresp2] = sift_corresp('img3.pgm','img2.pgm');

figure(2)
subplot(1,2,1)
%imshow(i3)
hold on;
plot(corresp1(:,2)',corresp1(:,1)','.');
subplot(1,2,2)
imshow(i2)
hold on;
plot(corresp2(:,2)',corresp2(:,1)','.');

no_of_iterations = 300;
dist_thresh = sqrt(10);
%Swapping Coords, so first column is x and second is y
% corresp1 = [corresp1(:,2)-size(i3,2)/2 corresp1(:,1)-size(i3,1)/2]; 
% corresp2 = [corresp2(:,2)-size(i2,2)/2 corresp2(:,1)-size(i2,2)/2];
corresp1 = [corresp1(:,2) corresp1(:,1)]; 
corresp2 = [corresp2(:,2) corresp2(:,1)];
for i=1:no_of_iterations
    ind = randsample(size(corresp1,1),4);
    A = [corresp2(ind(1),1) corresp2(ind(1),2) 1 0 0 0 -corresp1(ind(1),1)*corresp2(ind(1),1) -corresp1(ind(1),1)*corresp2(ind(1),2) -corresp1(ind(1),1);
         0 0 0 corresp2(ind(1),1) corresp2(ind(1),2) 1 -corresp1(ind(1),2)*corresp2(ind(1),1) -corresp2(ind(1),2)*corresp1(ind(1),2) -corresp1(ind(1),2);
         corresp2(ind(2),1) corresp2(ind(2),2) 1 0 0 0 -corresp1(ind(2),1)*corresp2(ind(2),1) -corresp1(ind(2),1)*corresp2(ind(2),2) -corresp1(ind(2),1);
         0 0 0 corresp2(ind(2),1) corresp2(ind(2),2) 1 -corresp1(ind(2),2)*corresp2(ind(2),1) -corresp2(ind(2),2)*corresp1(ind(2),2) -corresp1(ind(2),2);
         corresp2(ind(3),1) corresp2(ind(3),2) 1 0 0 0 -corresp1(ind(3),1)*corresp2(ind(3),1) -corresp1(ind(3),1)*corresp2(ind(3),2) -corresp1(ind(3),1);
         0 0 0 corresp2(ind(3),1) corresp2(ind(3),2) 1 -corresp1(ind(3),2)*corresp2(ind(3),1) -corresp2(ind(3),2)*corresp1(ind(3),2) -corresp1(ind(3),2);
         corresp2(ind(4),1) corresp2(ind(4),2) 1 0 0 0 -corresp1(ind(4),1)*corresp2(ind(4),1) -corresp1(ind(4),1)*corresp2(ind(4),2) -corresp1(ind(4),1);
         0 0 0 corresp2(ind(4),1) corresp2(ind(4),2) 1 -corresp1(ind(4),2)*corresp2(ind(4),1) -corresp2(ind(4),2)*corresp1(ind(4),2) -corresp1(ind(4),2)];
     H = reshape(null(A),[3,3])';
     [corresp1it,corresp2it] = deal(corresp1,corresp2);
     corresp1it(ind,:) = [];
     corresp2it(ind,:) = [];
     [x_rem,y_rem] = deal(corresp2it(:,1),corresp2it(:,2));
     [x_ac,y_ac] = deal(corresp1it(:,1)',corresp1it(:,2)');
     calculated = H*[x_rem';y_rem';ones(1,size(x_rem',2))];
     [x_calc,y_calc] = deal(calculated(1,:)./calculated(3,:),calculated(2,:)./calculated(3,:));
     dist = sqrt((x_calc-x_ac).^2+(y_calc-y_ac).^2);
     no_of_inliers = sum(dist<dist_thresh);
     if(no_of_inliers > 0.8*size(x_rem,1))
        break;
     end
end
figure(4)
%imshow(i3);
hold on
% plot(x_ac+size(i1,2)/2,y_ac+size(i1,1)/2,'.');
% hold on
% plot(x_calc+size(i1,2)/2,y_calc+size(i1,1)/2,'*');
plot(x_ac,y_ac,'.');
hold on
plot(x_calc,y_calc,'*');

H23 = H;
%% 
%Mosaicing
canvas = zeros(500,640*2.5);
canvas1 = canvas;
canvas2 = canvas;
H22 = [1 0 -size(canvas,2)/2;
       0 1 0;
       0 0 1];
H23 = eye(3);
for ii=1:size(canvas,1)
    for jj=1:size(canvas,2)
        i = ii-70;
        j = jj-size(canvas,2)/2 + size(i2,2)/2;%-size(canvas,2)/2;
        
%         i_1 = i;% - size(i1,1)/2;
%         j_1 = j;% - size(i1,2)/2;
%         tmp = H21*[j_1;i_1;1];
%         ii1 = tmp(1)/tmp(3);%+ size(i1,2)/2;
%         jj1 = tmp(2)/tmp(3);%+ size(i1,1)/2;
        %i1=i1+400;
        %j1=j1+180;
        
        i_3 = i;%- size(i1,1)/2;
        j_3 = j;%- size(i1,2)/2;
        tmp = H23*[j_3;i_3;1];
        ii3 = tmp(1)/tmp(3);%+ size(i1,2)/2;
        jj3 = tmp(2)/tmp(3);%+ size(i1,1)/2;
        
        i_2 = i;% - size(i2,1)/2;
        j_2 = j;% - size(i2,2)/2;
        tmp = eye(3)*[j_2;i_2;1];
        ii2 = tmp(1)/tmp(3);% + size(i2,2)/2;
        jj2 = tmp(2)/tmp(3);% + size(i2,1)/2;
           
        i_1 = i;%- size(i1,1)/2;
        j_1 = j;%- size(i1,2)/2;
        tmp = H21*[j_1;i_1;1];
        ii1 = tmp(1)/tmp(3);%+ size(i1,2)/2;
        jj1 = tmp(2)/tmp(3);%+ size(i1,1)/2;
        
        v2 = billinearH(ii2,jj2,[1 0 -size(canvas,2)/2; 0 1 0; 0 0 1],size(canvas),i2);
        v1 = billinearH(ii1,jj1,[1 0 -size(canvas,2)/2; 0 1 0; 0 0 1],size(canvas),i1);
        v3 = 0;
        %v3 = billinearH(ii3,jj3,[1 0 -size(canvas,2)/2; 0 1 0; 0 0 1],size(canvas),i3);
        %v3 = lab2bi(i2,H22,size(canvas));
        %v1 = billinearH(ii1,jj1,H21,size(canvas),i1);
        %v3 = billinearH(ii3,jj3,H23,size(canvas),i3);
        %if(jj==400)
        %    display(ii);
        %end
        %v2 = billinearH(i2,I,size(canvas));
        if((v1*v2>0)+(v1*v3>0)+(v2*v3>0) == 0)
            canvas(ii,jj) = int16(v1)+int16(v2)+int16(v3);
            canvas1(ii,jj) = 0;
        elseif((v1*v2>0)+(v1*v3>0)+(v2*v3>0) == 3) 
            canvas(ii,jj) = (int16(v1)+int16(v2)+int16(v3))/3;%v1+v3; %Overlap of 3
%             canvas(ii,jj)
%             (v1+v2+v3)/3
%             int16(v1)+int16(v2)+int16(v3)
%             v1
%             v2
%             v3
%             pause
            canvas1(ii,jj) = 128;
        else
            canvas(ii,jj) = (int16(v1)+int16(v2)+int16(v3))/2;%v1+v3; %Overalp of 2
            canvas1(ii,jj) = 255;
        end
        %canvas1(ii,jj) = v1;
        %canvas2(ii,jj) = v2;
    end
end
figure(5)
imshow(canvas,[0,255]);
figure(6)
imshow(canvas1,[0,255]);

        