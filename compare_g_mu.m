
% Load example data, which contains the following arrays
% t  : [nt, 1] vector of time values in seconds 
% x  : [nx, 1] vector of lateral positions meters
% v  : [nx, nt] array giving the vertical surface vibration velocity, as measured with OCE
% dt : [scalar] temporal sampling interval in seconds
% dx : [scalar] spatial sampling interval in meters 
% h  : [scalar] corneal thickness (in meters) used in the analytical forward model
load_data = importdata('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Cornea_iso_vacuum_Geqmu_fwP.txt');
data = load_data.data;

%load_ex_data = load('example_NITI_data.mat');
%test_data = load_ex_data.data_holder;



f_duplicates = data(:,1);
k_duplicates = data(:,2);
Vdb_flat = data(:,3);

f = unique(f_duplicates);
k = unique(k_duplicates);
VdB = reshape(Vdb_flat, length(f), length(k));
VdB = flipud(VdB');

f = f*1e3;

k = k*1e3;

fmax = 5000; % Hz
kmax = 2000; % 1/m
fmask = (f <= fmax);
kmask = (k <= kmax);
f = f(fmask);
k = k(kmask);
VdB = VdB(kmask, fmask);

%Input Parameters
h = 710e-6;
G0 = 26E3;
mu0 = 26E3;
reduce_flag = false;

G0_start = 1E3;
G0_end = 50E3;

mu0_start = 1E3;
mu0_end = 50E3;

num_points = 10;

G0_set = linspace(G0_start, G0_end, num_points);
mu0_set = linspace(mu0_start, mu0_end, num_points);

results = zeros(num_points, num_points,3);

if reduce_flag
    select_freq = [600, 900, 1200, 1500];
    num_freq = length(select_freq);
    
    % Initialize f_reduced as an array of selected frequencies
    f_reg = f;
    f_reduced = zeros(1, num_freq);
    
    % Loop over each frequency in select_freq
    for i = 1:num_freq
        [~, idx] = min(abs(f - select_freq(i)));  
        f_reduced(i) = idx;  
    end
    
    disp(f(f_reduced));
    
    f = f(f_reduced);
    VdB= VdB(:,f_reduced);
end
    
for n =1:num_points
    G0 = G0_set(n);
    for m =1:num_points
        mu0 = mu0_set(m);
        [G, mu, fval] = fit_spectrum_niti(f, k, VdB, h, G0, mu0);
        [mu_iso, fval_iso] = fit_spectrum_iso(f, k, VdB, h, G0);
        results(n,m,:) = [G, mu, mu_iso];
    end
end
[phi_max, ~] = global_mode_energy(f, k, VdB, 1);


gof_ti = abs(fval)/phi_max;
gof_iso = abs(fval_iso)/phi_max;


%
% Re-calculate NITI best fit dispersion curve
% The following assumptions are used inside the compute_niti_amode function.
% rho   : [scalar] assumed corneal tissue density 
% rho_l : [scalar] density of the liquid bounding the cornea (assume water)
% c_l   : [scalar] P-wave speed of the liquid bounding the cornea (assume water)
% cp    : [scalar] P-wave speed of the corneal tissue
% lambda : [scalar] stiffness tensor entry that enforces incompressibility, lambda >> mu > G
%
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
lambda = rho*cp^2 - 2*mu;
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
kfit = reshape(kfit, 1, []);


%
% Re-calculate isotropic best fit dispersion curve
% The following assumptions are used inside the fit_spectrum_niti function.
% rho   : [scalar] assumed corneal tissue density 
% rho_l : [scalar] density of the liquid bounding the cornea (assume water)
% c_l   : [scalar] P-wave speed of the liquid bounding the cornea (assume water)
% cp    : [scalar] P-wave speed of the corneal tissue
% lambda : [scalar] stiffness tensor entry that enforces incompressibility, lambda >> mu > G
%
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
lambda = rho*cp^2 - 2*mu_iso;
cfit_iso = compute_niti_amode(f, k, h, mu_iso, mu_iso, lambda, rho, rho_l, c_l);
kfit_iso = f(:)./cfit_iso(:);
kfit_iso = reshape(kfit_iso, 1, []);


%
% Plot results
%

save('kfit_data.mat','kfit');
save('kfit_iso_data.mat','kfit_iso');
%{
load_kfit = load('kfit_data.mat');
k_fit = load_kfit.kfit;
load_kfit_iso = load('kfit_iso_data.mat');
k_fit_iso = load_kfit_iso.kfit_iso;
%}

figure
pcolor(f/1000,k/1000, (VdB));
shading flat
%caxis([-70, -45])
cb = colorbar;
ylabel(cb, 'Power Spectrum (dB)')
xlabel('Frequency (kHz)')
ylabel('Wavenumber (1/mm)')
hold on
h_ti = plot(f(:)/1000, kfit(:)/1000, '-r', 'LineWidth', 1);
h_iso = plot(f(:)/1000, kfit_iso(:)/1000, '-k', 'LineWidth', 1);
%ylim()
legend([h_ti, h_iso], 'NITI', 'Iso', 'Location', 'NorthWest')
legend('boxoff')
colormap;
