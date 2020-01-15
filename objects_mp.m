
%% Preliminaries
clear all

fps = 10;  % frames/second

run_time = 3600;  % total run time

mp_cam = 12;  % megapixels (sensor size)

num_mp = fps * run_time * mp_cam;  % total megapixels over run time

obs_samp = 100;  % number of objects that comprise a sample

obj_size = linspace(1, 1000);  % object sizes
% obj_size = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480];

spectra = (10^6).*(4./3*pi*(obj_size./2).^3).^-1;  % Size spectra

px_bod = [10, 20, 50, 100];  % define the desired pixels per body size

beta = 3;  % acceptable blur size

lambda = 0.532;  % wavelength of light

% a place holder for the legend
leg_str = '%d px/diamter';
leg = cell(length(px_bod), 1);

%% 
obs_mp = zeros(size(obj_size,2), length(px_bod));
samp_hour = zeros(size(obj_size, 2), length(px_bod));

figure;
hold on

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

    % objects per megapixel
    obs_mp(:, ii) = v_frame.*spectra / mp_cam;
    
    leg{ii} = sprintf(leg_str, px_bod(ii));
    
    scatter(samp_hour(:, ii), obs_mp(:, ii))
    
end

%% plot formating
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('Sample Rate (n_{samples}/hour)', 'FontSize', 14)
ylabel('Objects/MP', 'FontSize', 14)
legend(leg(:), 'Location', 'northwest', 'FontSize', 12)


