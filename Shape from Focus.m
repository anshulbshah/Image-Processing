%% Loading the data and viewing the image loop. Initializing Variables
load stack.mat
for i = 1:numframes
    imshow(eval(['frame' sprintf('%.3d',i)]),[0,255])
end

delta = 50.50
%% Sum Modified Laplacian for each pixel each frame
focus_vals = zeros(115,115,numframes);
lap_vert = [0 1 0;
        0 -2 0;
        0  1 0];
lap_hor = [0 0 0;
         1 -2 1;
         0 0 0];
 nbd = 2;
 nbd_kernel = ones(2*nbd+1)
for i = 1:numframes     %Looping over all the frames
    image_frame = eval(['frame' sprintf('%.3d',i)]);
    sml = abs(convolution_operation(image_frame,lap_vert)) + abs(convolution_operation(eval(['frame' sprintf('%.3d',i)]),lap_hor));
    focus_vals(:,:,i) = convolution_operation(sml,nbd_kernel);
end
%% Estimating d_bar using shape from focus
depth_map = zeros(115,115);
d_vals = 0:delta:(numframes-1)*delta;
for l = 4:115
    for m = 4:115
        depth_map(l,m) = gaussian_interp(d_vals,reshape(focus_vals(l,m,:),numframes,1));
    end
end

