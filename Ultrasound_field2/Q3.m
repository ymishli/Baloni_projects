clc; close all;clear all;

% Initialize Field
%addpath('Field_II_ver_3_22_windows');
field_init(0)

% Generate the transducer aperture for send and receive
f0          = 3e6;              % Transducer center frequency [Hz]
fs          = 100e6;            % Sampling frequency [Hz]
c           = 1540;             % Speed of sound [m/s]
lambda      = c/f0;             % Wavelength [m]
width       = 0.29/1000;        % Width of element
height      = 5/1000;           % Height of element [m]
kerf        = 0.050/1000;       % Kerf [m]
focus       = [0 0 60]/1000;    % Fixed focal point [m]
N_elements  = 128;              % Number of physical elements
N_active    = 48;               % Active element on each
N_sub_x     = 1;                % Number of sub-divisions in x-direction of elements
N_sub_y     = 1;                % Number of sub-divisions in y-direction of elements
no_lines    = (N_elements-N_active)/2;               % Number of A-lines in image

% Set simulation parameters
set_sampling(fs);               % Sets sampling frequency
set_field('use_triangles',0);   % Tells whether to use triangles (1) or not (0)
set_field('use_rectangles',1);  % Tells whether to use rectangles (1) or not (0)
set_field('use_att',0);         % Tells whether to use attenuation (1) or not (0)
% set_field('c',c);             % Sets the speed of sound

% Generate aperture for transmission
tx = xdc_linear_array (N_elements, width, height, kerf, N_sub_x, N_sub_y, focus);
% Generate aperture for receive
rx = xdc_linear_array (N_elements, width, height, kerf, N_sub_x, N_sub_y, focus);
% Set apodization to hammimg widnow for tx and rx
guasian_apo = hamming(N_elements);
guasian_apo_matrix = repmat(guasian_apo' ,N_elements,1);
xdc_apodization(tx,zeros(N_elements,1),guasian_apo_matrix);
xdc_apodization(rx,zeros(N_elements,1),guasian_apo_matrix);
% Set the excitation of the transmit aperture
t = (0:1/fs:1.5/f0);
excitaion = sin(2*pi*f0*t);
xdc_excitation(tx,excitaion());
xdc_excitation(rx,excitaion());
% Set the impulse response
Bw = 0.6;
t_h = (-2/f0:1/fs:2/f0);
impulse_response = gauspuls(t_h,f0,Bw);    
xdc_impulse(tx,impulse_response);    
xdc_impulse(rx,impulse_response);    
% Set the focus point
xdc_focus(tx, 0, [0 0 40]/1000);
xdc_focus(rx, 0, [0 0 40]/1000);
% Do the calculation from one scatter point at focus
[v,t] = calc_scat_all(tx,rx,[0 0 40]/1000,1,1);
% Plot the individual responses
[N,M] = size(v);
scale = max(max(v));
v = v/scale;
figure('Name','Q3 - Section 1 - Individual traces');
timeline = (0:N-1)/fs+t;
for i = 1:128
    plot(timeline, v(:,i)+i,'b');
    hold on;
end
title('Q3 - Section 1 - Individual traces');
xlabel('time [sec]'); ylabel('Normalized response');
xlim([min(timeline) max(timeline)]);
ylim([1 128]);
hold off;
figure('Name','Q3 - Section 1 - Summed response');
plot((0:N-1)/fs+t, sum(v'));
title('Q3 - Section 1 - Summed response');
xlabel('time [sec]'); ylabel('Normalized response');
xlim([min(timeline) max(timeline)]);

% Do linear array imaging
dx=width; % Increment for image
z_focus=40/1000;
% Pre-allocate some storage
rf_data=zeros(1,no_lines);
for i=1:no_lines
% Find position for imaging
x=(i-1-no_lines/2)*dx;
% Set the focus for this direction
xdc_center_focus (tx, [x 0 0]);
xdc_focus (tx, 0, [x 0 z_focus]);
xdc_center_focus (rx, [x 0 0]);
xdc_focus (rx, 0, [x 0 z_focus]);
% Set the active elements using the apodization
apo=[zeros(1, 2*(i-1)) hamming(N_active)' zeros(1, N_elements-N_active-2*(i-1))];
xdc_apodization (tx, 0, apo);
xdc_apodization (rx, 0, apo);
% Calculate the received response
[v, t1]=calc_scat(tx, rx, [0 0 40]/1000, 1);
% Store the result
rf_data(1:max(size(v)),i)=v;
times(i) = t1;
end

[N,M] = size(rf_data);
rf_data = rf_data/max(max(rf_data));
figure('Name','Q3 - Section 2 - RF Data matrix result');
for i=1:no_lines  
    timeline = (0:N-1)/fs;    
    plot(timeline, rf_data(:,i)+i,'b');    
    hold on;
end
title('Q3 - Section 2 - RF Data matrix result');
xlabel('time [sec]'); ylabel('RF Data per Element');
hold off;

% build matrix aligned in time
times_shift = round((times - min(times))*fs);
figure('Name','Q3 - Section 3 - RF Data matrix aligned in time');
timeline = (0:N-1)/fs;    
for i=1:no_lines    
    rf_data(:,i) = circshift(rf_data(:,i),times_shift(i));
    plot(timeline, rf_data(:,i)+i,'b');    
    hold on;
end
title('Q3 - Section 3 - RF Data matrix aligned in time');
xlabel('time [sec]'); ylabel('RF Data per Element');
hold off;
figure('Name','Q3 - Section 3 - Received image');
min_sample=min(times)*fs;
depth=((0:size(N,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
imagesc(x*1000, depth*1000, rf_data);
colormap(gray(256));
title('Q3 - Section 3 - Gray scale image');
xlabel('Liteal axes [mm]'); ylabel('Depth [mm]');

% decimate by factor of 10
% rf_data_dec_10 = zeros(M/10,N);
dec = 10;
for i=1:no_lines    
    rf_data_dec_10(:,i) = decimate(rf_data(:,i),dec);        
end

% calc enelop with Hilbert transform
for i=1:no_lines
rf_env=abs(hilbert(rf_data(:,i)));
env(1:size(rf_env,1),i)=rf_env;
end

% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
figure('Name','Q3 - Section 4 -Image of point (40 dB dynamic range)');
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+40)/40;
depth=((0:size(env,1)-1)+min_sample)/(fs)*c/2;
x=((1:no_lines)-no_lines/2)*dx;
image(x*1000, depth*1000, env_gray);
xlabel('Lateral distance [mm]');
ylabel('Depth [mm]');
axis('image');
colormap(gray(128));
title('Q3 - Section 4 -Image of point (40 dB dynamic range)');

