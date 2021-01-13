clc; close all;

% Initialize Field
%addpath('C:\Users\Rony\Desktop\Ultrasound\Field2\Field_II_ver_3_22_windows');
field_init(0)

% Generate the transducer aperture for send and receive
f0          = 2.5e6;            %  Transducer center frequency [Hz]
fs          = 100e6;            %  Sampling frequency [Hz]
c           = 1490;             %  Speed of sound [m/s]
lambda      = c/f0;             %  Wavelength [m]
width       = 18.5/1000;        %  Width of element
height      = 13/1000;          % Height of element [m]
kerf        = 0;                % Kerf [m]
focus       = [0 0 60]/1000;    % Fixed focal point [m]
N_elements  = 1;                % Number of physical elements
N_sub_x     = 1;                % Number of sub-divisions in x-direction of elements
N_sub_y     = 1;                % Number of sub-divisions in y-direction of elements

% Set simulation parameters
set_sampling(fs);               % Sets sampling frequency
set_field('use_triangles',0);   % Tells whether to use triangles (1) or not (0)
set_field('use_rectangles',1);  % Tells whether to use rectangles (1) or not (0)
set_field('use_att',0);         % Tells whether to use attenuation (1) or not (0)
set_field('c',c);               % Sets the speed of sound

% Generate aperture for transmission
tx = xdc_linear_array (N_elements, width, height, kerf, N_sub_x, N_sub_y, focus);

% Set the excitation of the transmit aperture
t = (0:1/fs:1.5/f0);
excitaion = sin(2*pi*f0*t);
xdc_excitation(tx,excitaion);

%1.a
show_xdc_geir(tx,1);
axis equal;
view(3);
figure;
plot(t,excitaion);
title('excitation');
xlabel('t [sec]');
ylabel('V');

%1.b
% Set the impulse response of the transmit aperture
t_h = (-2/f0:1/fs:2/f0);
Bw = 0.6;
impulse_response = gauspuls(t_h,f0,Bw);
impulse_response = impulse_response.*sin(2*pi*f0*t_h);
xdc_impulse(tx,impulse_response);
figure;
plot(t_h,impulse_response);
figure;
freqz(impulse_response);

%1.c
% calculte the spatial impulse respponse of the transmit aperture
% [h_tx,start_time_h_tx] = calc_h(tx,[0,0,40]/1000);

%1.d
x = linspace(-10,10,100);
z = linspace(10,80,100);
y = zeros(1,100);
points = [x(:), y(:), z(:)];


