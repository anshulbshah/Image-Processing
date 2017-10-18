function [ d_peak ] = gaussian_interp( d_vals,curve )
    k = 3;
    F_mm1 = 0.0;
    F_m = 0.0;
    F_mp1 = 0.0;
    d_m = 0.0;
    delta = d_vals(2);
    while k<=size(d_vals,2)
        if(curve(k-1)>=F_m && curve(k-1)>=curve(k) && curve(k-1)>=curve(k-2))
            F_m = curve(k-1);
            F_mm1 = curve(k-2);
            F_mp1 = curve(k);
            d_m = d_vals(k-1);
        end
        k = k + 1;
    end
    d_mm1 = d_m-delta;
    d_mp1 = d_m+delta;
    %size(d_vals,2)
%     d_m
%     F_m
%     k
%     figure(1)
%     plot(d_vals,curve,'.');
%     hold on;
%     plot(d_m,F_m,'*');
%     hold on;
    d_peak = ((log(F_m)-log(F_mp1))*(d_m^2-d_mm1^2) - (log(F_m)-log(F_mm1))*(d_m^2-d_mp1^2))/(2*delta*((log(F_m)-log(F_mm1))+(log(F_m)-log(F_mp1))));
    %sigma_f = ((d_m^2-d_mm1^2)+(d_m^2-d_mp1^2))/(2*((log(F_m)-log(F_mm1))+(log(F_m)-log(F_mp1))));
    %peak = F_m/exp(-0.5*((d_m-d_peak)^2/sigma_f));
%     plot(d_peak,peak,'ro');
end

