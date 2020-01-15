% Model sample rate as a function of concentration and system mangnification
%
% Assumptions:
%   1. 12 MP sensor
%   2. Telecentric lens (implies orthographic projection)
%   3. Plankton only (no signal from detritus, etc)

clear all
close all

%% Preliminaries
fps = 10;  % frames per second
obs_samp = 100;  % number of objects that comprise a sample

obj_size = linspace(1, 1000, 100000);  % object sizes

% Size spectra assuming isometric scaling of metabolic rate and size abundance distribution of plankton
spectra = (10^6).*(4./3*pi*(obj_size./2).^3).^-1;  

% define the desired pixels per body size
px_bod = [10, 20, 50, 100];

% acceptable blur size
beta = 3;

% wavelength of light
lambda = 0.532;

% a place holder for the legend
leg_str = '%d px/diamter';
leg = cell(length(px_bod)+1, 1);
leg{1} = 'spectra';

%% compute the sample rate as a function of object size, concentration, and desired number of pixels per object
samp_hour = zeros(size(obj_size, 2), length(px_bod));

for ii = 1:length(px_bod)
    
    % minimum resolution needed to get desired number of pixels per body
    % length
    min_res = obj_size ./ px_bod(ii);
    
    % necessecary pixel size to achieve resolution (per body size)
    px_size = min_res ./ 2;
    
    % necessary numerical aperature (per body size)
    NA = 0.61*lambda ./ min_res;
    
    % mask out unrealistic NA
    NA(NA > 0.9) = 0.9;
    
    % max collection angle (per body size)
    theta_max = asin(NA);
    
    % depth of field (per body size)
    dof = beta.*px_size./tan(theta_max)./2;
    
    % volume of water per frame (per body size)
    v_frame = 4000 * 3000 .* dof .* px_size.^2 ./ (10000^3);
    
    % time to collect a sample (defined by SNR)
    t_per_samp = obs_samp ./ (v_frame .* spectra .* fps);
    
    % find the inflection point where no more objects are collected
    df1 = gradient(t_per_samp);
    df2 = gradient(df1);
    [~, ind] = max(df2);
    t_per_samp(1:ind) = 0;
    
    % compute the number of samples collected per hour
    samp_hour(:, ii) = 3600 ./ t_per_samp;
    
    leg{ii+1} = sprintf(leg_str, px_bod(ii));
    
end

%% plot it
figure;
yyaxis left
loglog(obj_size, spectra, 'LineWidth', 2);
ylabel('Concentration (n_obj/mL)', 'FontSize', 14)
yyaxis right
loglog(obj_size, samp_hour, 'LineWidth', 2);
ylabel('Sample rate (sample/hr)', 'FontSize', 14)
legend(leg(:), 'Location', 'southwest')
xlabel('Object size (um)', 'FontSize', 14)
