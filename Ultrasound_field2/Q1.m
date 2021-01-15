clc; close all;clear all;

% Initialize Field
%addpath('C:\Users\Rony\Desktop\Ultrasound\Field2\Field_II_ver_3_22_windows');
field_init(0)

% Control the plots plot_on = [sec_a sec_b sec_c sec_d sec_e]
plot_on = [0 0 1 0 0];

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
xdc_excitation(tx,excitaion());

%1.a
if plot_on(1) == 1    
    figure('Name','Q1 - Section 1 - Transducer Geometry');
    show_xdc_geir(tx,1);    
    axis equal;
    view(3);
    title('Q1 - Section 1 - Transducer Geometry');
    figure('Name','Q1 - Section 1 - Excitation signal');
    plot(t,excitaion);
    title('Q1 - Section 1 - Excitation signal');
    xlabel('t [sec]');
    ylabel('V');
end

%1.b
Bw = 0.6;
t_h = (-2/f0:1/fs:2/f0);
impulse_response = gauspuls(t_h,f0,Bw);    
xdc_impulse(tx,impulse_response);    
if plot_on(2) == 1    
    figure('Name','Q1 - Section 2 - Impulse Response in Time');
    plot(t_h,impulse_response);    
    title('Q1 - Section 2 - Impulse Response in Time');
    xlabel('t [sec]');
    ylabel('V');
    figure('Name','Q1 - Section 2 - Impulse Response in Frequncy domain');    
    freqz(impulse_response);
    title('Q1 - Section 2 - Impulse Response in Frequncy domain');
end

%1.c
if plot_on(3) == 1
    [h, start_time] = calc_hp (tx,[0,0,40/1000]);
    figure;
    plot(h);
end


%1.d
if plot_on(4) == 1
    [x,y,z]=meshgrid(linspace(-10,10,100)/1000,0,linspace(10,80,100)/1000);
    points=[x(:) y(:) z(:)];
    Im_size=[length(x),1,length(z)];
    [hp,start_t]=calc_hp(tx,points);
    [m,n]=size(hp);

    % With 'Norm' on each impulse response
    for i=1:n
      P1(i) = norm(hp(:,i));
    end
    P1=reshape(P1,[Im_size(1),Im_size(3)]);
    P1=rot90(P1,1);
    Result=flipud(P1);
    Result_norm=Result-min(min(Result));
    Result_norm=Result_norm/max(max(Result_norm));
    figure;
    subplot(1,2,1);
    imagesc(1000*x(:),1000*z(:),Result_norm);
    colormap(hot)
    title('Transmit Field Example');
    xlabel('X[mm]');ylabel('Z[mm]');
    % Convert Intensity results to dB
    %  Result_dB=Convert2dB(Result);
    %  subplot(1,2,2);
    %  imagesc(1000*x(:),1000*z(:),Result_dB);
    %  colormap(hot);
    %  title('Transmit Field Example dB');
    %  xlabel('X[mm]');ylabel('Z[mm]');
    %  colorbar
end
