function interactive_spatial_plot(data,slice,samplesPerSec)
     
    
    
    if nargin <1 
    %Suture
    %DataName = '20210920_TD_SutureTriP_20ns_20ns_100res_200GV_200kHz_200Avg_1f_3fpm.mat';
    %DataPath = '/Users/calvinsmith/Bouma_lab/PAT_data/TD_SutureTriP_20ns_20ns_100res_200GV_200kHz_200Avg_1f_3fpm/PA/';

    %Adipocyte
    %DataName = '20210730_TD_BAT_dish2_ROI11c_1192_100res_30GV_100kHz_500Avg_1f_1fpm.mat';
    %DataPath = '/Users/calvinsmith/Bouma_lab/PAT_data/TD_BAT_dish2_ROI11c_1192_100res_30GV_100kHz_500Avg_1f_1fpm/PA/';

    %Adipoctye 2
    %DataName = '20210728_TD_WAT_dish2_441res_200x200FOV_100kHz_200Avg.mat';
    %DataPath = '/Users/calvinsmith/Bouma_lab/PAT_data/WAT_dish2/TD_WAT_dish2_441res_200x200FOV_100kHz_200Avg/PA/';
    
    %Adipocyte 3
    DataName = '20210728_TD_WAT_dish2_Align_100res_30GV_100kHz_50Avg_1f_5fpm.mat';
    DataPath = '/Users/calvinsmith/Bouma_lab/PAT_data/WAT_dish2/TD_WAT_dish2_Align_100res_30GV_100kHz_50Avg_1f_5fpm/PA/';
    
    %Custom Data
   
    raw_data = load([DataPath, DataName]);
    data = raw_data.S;

    
    end
    if nargin <2 
    slice = 150;
    end
    if nargin < 3
    samplesPerSec = raw_data.samplesPerSec;
    end


    cleanupObj = onCleanup(@() close_all_figures);

    function close_all_figures
        close all;
        delete(findall(0, 'Type', 'figure'));     % Closes all `figure` windows
        delete(findall(0, 'Type', 'uifigure'));
    end

    nF = 1;
   
    f_BP = [5 45]*1e6;
   
    

    chebyshev_data = zeros(size(data));
    

    Nx = sqrt(size(data,1));
    Ny = Nx;
    NPoints = Nx*Ny;
    N_time_samples = size(data,2);


    chebyshev_data = chebyshev_filter(data,Nx,Ny,NPoints,N_time_samples);
    chebyshev_data = construct_bidirectional_scan(chebyshev_data,Nx,Ny);
    chebyshev_data = reshape(chebyshev_data,Nx,Ny,N_time_samples,nF);

    filtered_data = full_filter_data(data,NPoints,f_BP,nF,samplesPerSec);
    filtered_data = construct_bidirectional_scan(filtered_data,Nx,Ny);
    filtered_data = reshape(filtered_data,Nx,Ny,N_time_samples,nF);

    hilbert_data = abs(hilbert(filtered_data));
    %hilbert_data = denoise_data(hilbert_data);
    
    custom_data = imgaussfilt3(hilbert_data,2);

    data = construct_bidirectional_scan(data,Nx,Ny);
    data = reshape(data,Nx,Ny,N_time_samples,nF);
   

    magnitude_data = zeros(size(data));
    magnitude_data_raw = get_intensity_magnitude_data(filtered_data,N_time_samples);

    magnitude_data(:,:,150,1) = magnitude_data_raw;
    
    Nx = size(data,1);
    Ny = size(data,2);
    N_time_samples = size(data,3);
    NPoints = Nx*Ny;


    button_f = figure;
    current_data = data;
    open_flag = true;

        f = figure('Position',[100,100,1200,400]);
        zoom on;
        pan on;
        setappdata(f, 'xCoord', 1); % Initial default X value
        setappdata(f, 'yCoord', 1); % Initial default Y value
        ax1 = subplot(1, 3, 1);
       
        img = imagesc([], [], data(:,:,slice,nF)); 
        
        colormap hot;
        colorbar;
        title(ax1, 'Original Data');
        zoom on;
        pan on;
        menu_x_cord = 100;
        menu_y_cord = 300;
        spacing = 40;
        
        guidata(f, current_data);  

        cheby_toggle_btn = uicontrol( button_f, 'Style', 'togglebutton', 'String', 'Switch to Chebyshev', ...
                                'Position', [menu_x_cord, menu_y_cord , 150, 30], 'Value', 0, ...
                                'Callback', @(src, event)  toggleData(src, img, data, chebyshev_data, slice, nF, 'Chebyshev',f,ax1));
        filter_toggle_btn = uicontrol( button_f, 'Style', 'togglebutton', 'String', 'Switch to Filtered', ...
                                'Position', [menu_x_cord, menu_y_cord - spacing, 150, 30], 'Value', 0, ...
                                'Callback', @(src, event)   toggleData(src, img, data, filtered_data, slice, nF, 'Filtered',f,ax1));
        hilbert_toggle_btn = uicontrol( button_f, 'Style', 'togglebutton', 'String', 'Switch to Hilbert', ...
                                'Position', [menu_x_cord, menu_y_cord - 2*spacing, 150, 30], 'Value', 0, ...
                                'Callback', @(src, event)  toggleData(src, img, data, hilbert_data, slice, nF, 'Hilbert',f,ax1));
        magnitude_toggle_btn = uicontrol( button_f, 'Style', 'togglebutton', 'String', 'Switch to Magnitude(Final)', ...
                                'Position', [menu_x_cord, menu_y_cord - 3*spacing, 150, 30], 'Value', 0, ...
                                'Callback', @(src, event)  toggleData(src, img, data, magnitude_data, slice, nF, 'Magnitude(Final)',f,ax1));
        
        custom_toggle_btn = uicontrol( button_f, 'Style', 'togglebutton', 'String', 'Switch to CUSTOM', ...
                                'Position', [menu_x_cord, menu_y_cord - 4*spacing, 150, 30], 'Value', 0, ...
                                'Callback', @(src, event)  toggleData(src, img, data, custom_data, slice, nF, 'CUSTOM',f,ax1));
       
        plot_vol_btn = uicontrol( button_f, 'Style', 'pushbutton', 'String', 'PlotVolumeData', ...
                                'Position', [menu_x_cord, menu_y_cord-5*spacing , 150, 30], ...
                                        'Callback', @(src, event) plot_vol_data(guidata(f)));
                
        uicontrol(button_f, 'Style', 'text', 'String', 'X Coord:', ...
                  'Position', [menu_x_cord + 200, menu_y_cord, 60, 20]);
        x_input = uicontrol(button_f, 'Style', 'edit', 'Position', [menu_x_cord + 260, menu_y_cord, 60, 20]);
        
        
        uicontrol(button_f, 'Style', 'text', 'String', 'Y Coord:', ...
                  'Position', [menu_x_cord + 200, menu_y_cord - spacing, 60, 20]);
        y_input = uicontrol(button_f, 'Style', 'edit', 'Position', [menu_x_cord + 260, menu_y_cord - spacing, 60, 20]);
        
        
        uicontrol(button_f, 'Style', 'pushbutton', 'String', 'Select Point', ...
          'Position', [menu_x_cord + 200, menu_y_cord - 2*spacing, 120, 30], ...
          'Callback', @(src, event) updateCoordinates(x_input, y_input, f));

        
      
        while true
            
            figure(f);
            x = 30;
            y= 40;
            [x,y,button]=ginput(1);



            %if button ~= 1
             %       break;
            %end
            subplot(1,3,2);
           
            cla;
            hold on
            x = round(x);
            y = round(y);


            current_data = guidata(f);
            data = current_data;
         
            %b = uibutton(point_f,"state");
            point = data(x,y,:,nF);
            point = squeeze(point);
            
            plot(point)
            
            
            plot_hilbert = plot(point, color = 'red');
            plot_hilbert.Visible = 'off';
           
        
            title_text = sprintf("Spatial Point %d %d vs time",x,y);
            title(title_text)
            legend({'Spatial Data, Hilbert Transformed Data'});
            
           
            subplot(1,3,3)
            hold on;
            cla;
           
            fft_point = abs(fftshift(fft(point)));
            plot_fft = plot(fft_point);
            ylim([0 500]);
            
             
            corr_point = xcorr(point);
            plot_corr = plot(corr_point);
            plot_corr.Visible = 'off';

            title("FFT and Others")
            btn = uicontrol('Parent',button_f, 'Style', 'togglebutton', 'String', 'Hilbert', ...
                'Position', [460 menu_y_cord 100 30], 'Value', 0, ...
                'Callback', @(src, event) toggleFilter(src, plot_hilbert));
            fft_btn = uicontrol('Parent',button_f, 'Style', 'togglebutton', 'String', 'FFT', ...
                'Position', [460 menu_y_cord - spacing 100 30], 'Value', 0, ...
                'Callback', @(src, event) toggleFilter(src, plot_fft));
            corr_btn = uicontrol('Parent',button_f, 'Style', 'togglebutton', 'String', 'Auto_corr', ...
                'Position', [460 menu_y_cord - 2*spacing 100 30], 'Value', 0, ...
                'Callback', @(src, event) toggleFilter(src, plot_corr) );
            

        end  
end

function toggleFilter(button, plot_hilbert)
    % Toggle the visibility of the Hilbert transform plot
    if button.Value == 1
        plot_hilbert.Visible = 'on';
    else
        plot_hilbert.Visible = 'off';
    end
end

function updateCoordinates(x_input, y_input, fig_handle)
    
    x_str = get(x_input, 'String');
    y_str = get(y_input, 'String');
 
    x = str2double(x_str);
    y = str2double(y_str);    
 
    x = round(x);
    y = round(y);
    setappdata(fig_handle, 'xCoord', x);
    setappdata(fig_handle, 'yCoord', y);

    
    figure(fig_handle); 
    drawnow;    
end

function [current_data] = toggleData(button, img, data, transform_data, slice, nF, transform_type,fig_handle,ax_handle)
    % Switch between original data and Chebyshev data on the heatmap
    if button.Value == 1
        img.CData = transform_data(:, :, slice, nF); 
        button.String = 'Switch to Original';
        current_data = transform_data;
        title(ax_handle, sprintf('%s Data', transform_type));
    else
        img.CData = data(:, :, slice, nF); 

        button.String = sprintf('Switch to %s',transform_type);
        current_data = data;
         title(ax_handle, 'Original Data');
    end
    guidata(fig_handle, current_data);
end

%%Helper Functions

function S = chebyshev_filter(S,Nx,Ny,NPoints,N_time_samples)
[b, a] = cheby1(4, .01, 2*0.01, 'high');                                         % Calculate highpass Chebyshev filter coeff. of order 4

NySplit = 1;

% Chunks data by NySplit
% CS: This is relevant is NySplit != 1
for i = 1:NySplit                                                                % Split filtering length to avoid 'Out of memory'
    % iStart is the index of the element which starts for the
    % NySplit chunk i
    % Same for iEnd just the index of the end element for the
    % NySplit chunk i
    iStart = Ny/NySplit*(i-1)*Nx + 1;
    iEnd   = min(Ny/NySplit*i*Nx, NPoints);

    for j = 1:N_time_samples                                                           % Loop over samples
        ss = squeeze(double(S(iStart:iEnd, j)));                                 % Transform matrix S(s:f,j) into vector [s f]
        S(iStart:iEnd, j) = single(filtfilt(b, a, ss));                          % Apply 0-phase Chebyshev filter of order 8 to samples j
    end
end
end

function S = full_filter_data(S,NPoints,f_BP,nF,samplesPerSec)

    chunk_size  = 4000;   
    N_time_samples = size(S,2);
     trigDelay = 1152;
    dt = 1/samplesPerSec;
        t = 0:dt:(N_time_samples-1)*dt;
        t = t + trigDelay*dt; 
    
    % Number of signals to be filtered simultaneously, chunk size for filtering
        % Number of chunks to be filtered
    NIterations = ceil(NPoints/chunk_size);
        f = linspace(-samplesPerSec/2, samplesPerSec/2, length(t)); f = f';

        %Creates a sixth-order Gaussian and replicates for N Rows to
        %filtArray which is used later a the bandbass filter

        filtOrder = 6;
        filtArray = (exp(-(f/f_BP(2)).^filtOrder)) ...
            .* (1-(exp(-(f/f_BP(1)).^filtOrder)));
        filtArray = repmat(filtArray, 1, chunk_size);

            
        S = permute(S,[2 1 3]);
        S = filter_data(S, chunk_size,NPoints,NIterations,nF,filtArray);
end


function S = fft_filter(S, samplesPerSec,f_BP)

%NaScan  = 10000;                                                                  % define the maximum number of measurements that should be filtered at once (to prevent excessive use of memory)
%CS: I'm going to change NaScan to 1000
NaScan = 1000;

% How many time it has to split S
nSplit  = ceil(size(S,1)/NaScan);
fprintf("nSplit")
disp(nSplit)
%Size(S,2) = N_time_samples, Size(S,1) = NPoints
%filtPat is an 2D array of ones for each measurement by each sample
filtPat = ones(NaScan, size(S,2), 'single');
%NNtime is an 1D array from 1 -> N_time_samples
NNtime  = 1:size(S,2);

% Filterwidth details
if samplesPerSec == 1e9                                                                     % the larger the filter width the broader are the stripes that are removed
    FiltWidth = 50;                                                              % Be aware, that large (low-frequency) structures are biased if the filter is too broad
elseif samplesPerSec == 500e6
    FiltWidth = 25;
else
    FiltWidth = 37.5;
end

% for high frequency reconstruction the FFT-filter can be set more
% restrictive (filters out broad stripes as well)
if f_BP(1) > 90e6
    FiltWidth = FiltWidth * 4;
elseif f_BP(1) > 30e6
    FiltWidth = FiltWidth * 2;
elseif f_BP(1) > 10e6
    FiltWidth = FiltWidth * 1.5;
end

% the first exponential expression is a Gaussian centered at the
% beginng of NNtime while the second is centered at the end of
% NNtime
% By summing them you create a "two-humped" Gaussian around the end
% and beginning which you could use to filter out low/high
% frequencies
filtPat(1,:) = exp(-NNtime.^2./(2*(length(NNtime)/FiltWidth)^2))...              % filter out vertical stripes in b-scan by removing low frequency components of the central slice (along y) in the FFT image
    + exp(-(NNtime-length(NNtime)).^2./(2*(length(NNtime)/FiltWidth)^2));

for iSplit = 1:nSplit
    if iSplit == nSplit
        % Grabs the remaining points
        if nSplit ~= 1
            nRest        = size(S((nSplit-1)*NaScan:end,:));
        else
            nRest        = size(S((1)*NaScan:end,:));
            nSplit = 2;
        end
        % Changes filtPat size to nRest instead of N_time_samples
        filtPat      = ones(nRest,'single');
        filtPat(1,:) = exp(-NNtime.^2./(2*(length(NNtime)/FiltWidth)^2))...
            + exp(-(NNtime-length(NNtime)).^2./(2*(length(NNtime)/FiltWidth)^2));
        % Filters frequens using filtPat after FFT and then uses
        % inverse FFT to return the filtered image
        fftImage     = fft2(S((nSplit-1)*NaScan:end,:)).*filtPat;
        ifftImage    = ifft2(fftImage);

        S((nSplit-1)*NaScan:end,:) = real(ifftImage);
    else
        %Filters frequencies using filtPat for all the samples in
        %each measurement foreach point in S
        % CS: Need to check this reasoning
        fftImage     = fft2(S((iSplit-1)*NaScan+[1:NaScan],:)).*filtPat;
        ifftImage    = ifft2(fftImage);

        % CS: Should probably create a new matrix instead of
        % reassgning elements to S
        % Reassigns elements of S to the real components of the
        % produced image
        S((iSplit-1)*NaScan+[1:NaScan],:) = real(ifftImage);
    end
end

clear filtPat fftImage ifftImage;
end

function cut1 = plot_data_get_cut(S)

figure;
imagesc([], [], S(:,:,1)); colormap gray;
[cut1,cut2]=ginput;
disp(cut1)
cut1 = int16(cut1);
end

function S = hilbert_transform(S,NPoints)
S = double(S);
for i=1:NPoints
    % For each point i performs the Hilbert Transfrom on the data at i
    % and uses absolute value to calculate the envelope of the
    % frequency
    % It stores this new envelope data in S(i, :, :)
    S(i,:,:) = abs(hilbert(S(i,:,:)));
end
end

function S = construct_bidirectional_scan(S,Nx,Ny)

for i = 1:Ny
    if rem(i,2) == 0
        s = S((i-1)*Nx+1:i*Nx,:,:);

        S((i-1)*Nx+1:i*Nx,:,:) = s(end:-1:1,:,:);
    end
end
end


function result = plot_image_hot(H,rotFlag,c,ff)
figure;
%Plots H and rotated by rotFlag
plot_data = H(:,:,c,ff);
plot_data = rot90(plot_data,rotFlag);

%imagesc([Y(end)/1e3 Y(1)/1e3 ], [X(end,end)/1e3 X(1,1)/1e3], plot_data);
%imagesc([Y(end) Y(1)]/1e3, [X(end) X(1,1)]/1e3, plot_data);
imagesc([],[], plot_data);
colormap hot;
colorbar;

%        colorbar('southoutside');
xlabel(gca, 'y (mm)', 'FontWeight', 'bold', 'units', 'normalized', 'FontSize', 12);
ylabel(gca, 'x (mm)', 'FontWeight', 'bold', 'units', 'normalized', 'FontSize', 12);
set(gca,'YDir','normal');
axis image
%         set(gcf, 'Position', get(0, 'Screensize'));

end


function result = plot_image(H,c,title_text)
if nargin < 3
    title_text = "AAA";
end
ff = 1;
rotFlag = 0;
figure;

plot_data = H(:,:,c,ff);
plot_data = rot90(plot_data,rotFlag);

%imagesc([Y(end)/1e3 Y(1)/1e3 ], [X(end,end)/1e3 X(1,1)/1e3], plot_data);
%imagesc([Y(end) Y(1)]/1e3, [X(end) X(1,1)]/1e3, plot_data);
imagesc([],[], plot_data);
title(title_text);
colormap hot;
colorbar;
end

function S = filter_data(S, chunk_size,NPoints,NIterations, nF,filtArray)
for ii = 1:NIterations
            %Creates an array SS which selects a chunk of chunk_size from S for
            %each iteration ii
            SS = S(:, (ii-1)*chunk_size+1:min(ii*chunk_size, NPoints),:);
            % CS: Confused on how this works
            SS = fftshift(fft(ifftshift(SS, 1), [], 1), 1);                          % Take fourier transform of the signals (FFT)

            if NIterations == ii
                filtArray = filtArray(:, 1:size(SS, 2));
            end

            % This applies the bandpass filter in the frequency domain?
            % CS: Need to think about this abit more?
            for ff = 1:nF;
                SS(:,:,ff) = real(fftshift(ifft(ifftshift(SS(:,:,ff).*filtArray, 1), [], 1), 1));        % Multiply with f-window and transform back to time domain (iFFT)
            end
            % CS: Again weird reassigning to S instead of making a new data
            % structure
            % Reassigns SS to the "S Chunk of chunk_size " with the filtered
            % data
            S(:, (ii-1)*chunk_size+1:min(ii*chunk_size, NPoints),:) = SS;

        end
        % Rearranges S so it is back to  S(NPoints, N_time_samples, aAvg?)
        S = permute(S,[2 1 3]);

        clear filtArray SS;
        disp('Done Filter')
end



function result = plot_image_bidirect(S,Nx,Ny,N_time_samples,nF,title_text)
if nargin < 6
    title_text = "";
end
S_test = construct_bidirectional_scan(S,Nx,Ny);
S_test = reshape(S_test,Nx,Ny,N_time_samples,nF);
plot_image(S_test,100,title_text)
end


function result = plot_spatial_point_v_time_bidirect(x,y,S,Nx,Ny,N_time_samples,nF,title_text)
if nargin < 6
    title_text = "";
end
figure;
S_test = construct_bidirectional_scan(S,Nx,Ny);
S_test = reshape(S_test,Nx,Ny,N_time_samples,nF);
point = S_test(x,y,:);
point = squeeze(point);
plot(point)
title_text = sprintf("Spatial Point %d %d vs time",x,y);
title(title_text)
end


function result = plot_spatial_point_v_time(x,y,H,nF,title_text)
if nargin < 6
    title_text = "";
end
figure;
point = H(x,y,:,nF);
point = squeeze(point);
plot(point)
title_text = sprintf("Spatial Point %d %d vs time",x,y);
title(title_text)
end


function H = get_intensity_magnitude_data(R_data,N_time_samples)
            
cut1 = [1 N_time_samples];
c = 1;
ff = 1;
plotAxis = 3;
R_selected_data = R_data(:,:,cut1(c):cut1(c+1),ff);
R_max_xy_plane = double(squeeze(max(R_selected_data, [], plotAxis)));
R_min_xy_plane = double(squeeze(min(R_selected_data, [], plotAxis)));

H(:,:,c,ff) = R_max_xy_plane - R_min_xy_plane;
end




function plot_vol_data(R_data)
    if nargin < 1
                
        DataName = 'r_TD_SutureTriP_20ns_20ns_100res_200GV_200kHz_200Avg_1f_3fpm.mat';
        DataPath = '/Users/calvinsmith/Bouma_lab/PAT_data/TD_SutureTriP_20ns_20ns_100res_200GV_200kHz_200Avg_1f_3fpm/PA/';
        
        data = load([DataPath, DataName] );
        R_data = data.R_data;

    end

    
    H_vol = squeeze(R_data(:,:,:,1));
    H_vol_max = max(H_vol,[],'all');
    H_vol_min = min(H_vol,[],'all');
    disp(H_vol_max);
    disp(H_vol_min)

    threshold = prctile(H_vol(:), 96);  %

    
    disp([' intensity threshold: ', num2str(threshold)]);

   
    intensity = linspace(H_vol_min, H_vol_max, 256);
    
    
    alpha_map = zeros(size(intensity));
    alpha_map(intensity > threshold) = 0.8;  % Full opacity for top 30%
    
    %color = [0 0 0; 0.5 0.5 0.5; 1 1 1];  % Example grayscale color map
    %color_map = interp1([H_vol_min, threshold, H_vol_max], color, intensity);
    color_map = hot(256);

    sx = 1;
    sy= 1;
    sz = 1.3;
    A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
    tform = affinetform3d(A);
    volshow(H_vol, Alphamap = alpha_map, Colormap = color_map, Transformation=tform)
    
end

function clean_R_data = denoise_data(R_data)
    net = denoisingNetwork("DnCNN");
    R_data = permute(R_data, [1 2 4 3]);
    
    clean_R_data = denoiseImage(R_data,net);
    clean_R_data = permute(clean_R_data, [1 2 4 3]);
end




