function c_mask =  find_concentration_mask(c_data, hb_S, hbo_S)
    if nargin <1 
        c_data = load("total_concentration_recon.mat");
        c_data = c_data.sum_C
    end
    if nargin < 2  
        hb_S = load('hb_S_data.mat','hb_S');
        hb_S = hb_S.hb_S;
    end
    if nargin < 3
        hbo_S = load('hbo_S_data.mat','hbo_S');
        hbo_S = hbo_S.hbo_S;
    end




mask = load('c_mask.mat');
mask= mask.mask;

figure;
subplot(3,2,1)
imagesc(c_data);
title('Concentration Sum')
colormap;
colorbar;

subplot(3,2,2)
imagesc(hb_S);
title('HB Saturation')
colormap;
colorbar;

subplot(3,2,3)
imagesc(hbo_S);
title('HBO Saturation')
colormap;
colorbar;


gradient_mask = abs((gradient(c_data)));
subplot(3,2,4)
imagesc(gradient_mask);
title('Gradient Mask')
colormap;
colorbar;

subplot(3,2,5)
imagesc(mask)
title('Real Mask')
colormap;
colorbar;


masked_hbo_S = mask .* hbo_S;

subplot(3,2,6)
imagesc(masked_hbo_S)
title('HbO Saturation Mask')
colorbar;


% Weighted average

Nx = size(c_data,1);
Ny = size(c_data,2);

w_avg_hbo_S = zeros(size(c_data));

w_avg_hbo_S = hbo_S.*(c_data.^2);

w_avg_score =sum(w_avg_hbo_S)/(sum(c_data.^2));

figure;
subplot(1,2,1)
imagesc(w_avg_hbo_S)
title('weighted average');
colormap;
colorbar;

subplot(1,2,2)
imagesc(masked_hbo_S)
title('Masked Hbo S');
colormap;
colorbar;








%{
for i = 1:Nx
    for j = 1:Ny
        w_avg_hbo_S(i,j) = hbo_S(i,j)c_data(i,j)
    end
end
%}



end


