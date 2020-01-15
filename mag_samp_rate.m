% Model sample rate as a function of concentration and system mangnification
%
% Assumptions:
%   1. 12 MP sensor
%   2. Telecentric lens (implies orthographic projection)
%   3. Plankton only (no signal from detritus, etc)

clear all
close all

%% Microscope preliminaries
% Define the desired microscope magnifications
mag = [5, 2, 1, 0.5, 0.25, 0.138, 0.05];

% the associated pixel size in object space (assuming 12 MP sensor)
px = 3.45./mag;

% Numerical aperature (assuming 532 nm illumination and 2 px resolution)
NA = 0.61*(0.532./(2*px));

% maximum collection angle
theta_max = asin(NA);

% depth of focus (assuming a constant acceptable blur)
dof = (3*px)./(tan(theta_max));

% sample volume in uL (the amount of water imaged in an individual frame)
v_frame = (px.^2 .* dof * 4000 * 3000 * 1000) ./ 10000^3;

%% Size distribution and concentration
% Define the assumed concentration of plankton (in microns)
obj_size = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480];

% compute concentration per uL
con = 1000000 * obj_size.^-3;

%% Compute the resolution for each body size based on the magnification (number of pixels per body length)
res = zeros(size(obj_size, 2), size(px, 2));

for ii = 1:length(px)
    res(:, ii) = obj_size / px(ii);
end

%% Compute the sample rate
fps = 8;  % frames per sec

f_hour = fps * 3600;  % number of frames per hour

obs_hour = zeros(size(obj_size, 2), size(res, 2));

% compute the total number of objects per hour for each microscope
% objective
for ii = 1:length(px)
    obs_hour(:, ii) = con * f_hour * v_frame(ii) / 1000;
end

% compute the number of samples, where sample is defined by the desired SNR
samp = 100;  % the desired number of objects to get a statistically relevant collection of objects

samp_rate = obs_hour ./ samp;  % rate of sample collection

%% Mask the sample rate based on desired system parameters

min_size = 20;  % minimum number of pixels per object
min_samp = 0.001;  % minimum allowable samples per hour

size_mask = (res > min_size);  % mask based on system resolution

rate_mask = (samp_rate > min_samp);  % mask based on desired SNR

mask = size_mask .* rate_mask;

% mask the sample rate matrix
samp_rate = samp_rate .* mask;

%% plot it
leg = {'5x', '2x', '1x', '0.5x', '0.25x', '0.138x', '0.05x'};
figure;
loglog(obj_size, samp_rate, 'LineWidth', 2)
xlabel('Object Size (um)', 'FontSize', 14)
ylabel('Sample Rate (samp/hour)', 'FontSize', 14)
ylim([0.01, 200])
legend(leg(:))
