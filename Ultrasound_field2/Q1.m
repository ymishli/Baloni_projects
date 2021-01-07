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

