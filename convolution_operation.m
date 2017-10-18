function [ image_f ] = convolution_operation(i2,kernel)
%Kernel size has to be odd
    i2 = im2double(i2);
    size_o = floor(size(kernel,1)/2);
    image_f = zeros(size(i2));
    for l = (size_o+1):(size(image_f,1)-size_o)
        for m = (size_o+1):(size(image_f,2)-size_o)
            image_f(l,m) = sum(sum(times(kernel,i2(l-size_o:l+size_o,m-size_o:m+size_o))));
        end
    end
end

