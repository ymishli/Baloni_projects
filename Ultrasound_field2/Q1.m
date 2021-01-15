clc; close all;clear all;

% Initialize Field
%addpath('C:\Users\Rony\Desktop\Ultrasound\Field2\Field_II_ver_3_22_windows');
field_init(0)

% Control the plots plot_on = [sec_a sec_b sec_c sec_d sec_e]
plot_on = [0 0 0 1 0];

% Generate the transducer aperture for send and receive
% f0          = 2.5e6;            %  Transducer center frequency [Hz]
f0          = 5e6;            %  Transducer center frequency [Hz]
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
% set_field('c',c);               % Sets the speed of sound

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
    [h, start_time] = calc_hp(tx,[0,0,40/1000]);
    N = length(h);
    t_h_ext = 1/fs*(0:1:N-1)- t_h(1);
    figure('Name','Q1 - Section 3 - Field Reaction on point [0,0,40mm]');
    subplot(2,1,1);
    plot(t_h_ext,h, 'Color', 'g');
    title('Q1 - Section 3 - Field Reaction on point [0,0,40mm]');
    xlabel('t [sec]');
    ylabel('V');
    h_norm = h/max(h);
    impulse_response_norm = impulse_response/max(impulse_response);    
    subplot(2,1,2);
    p = plot(t_h,impulse_response_norm, t_h_ext+start_time,h_norm);
    p(1).Color = 'b';
    p(2).Color = 'g';    
    title('Q1 - Section 3 - Normalized Field Reaction on point [0,0,40mm] - relative to Excitation');
    annotation('textbox', [0.3, 0.2, 0.5, 0.1], 'String', "the time diffrence between excitation and his response = 40mm/c");
    xlabel('t [sec]');
    ylabel('V');
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
    figure('Name','Q1 - Section 4 - Transmit field picture at [-10,10]mm*[10,80]mm');
    subplot(1,2,1);
    imagesc(1000*x(:),1000*z(:),Result_norm);
    colormap(hot)
    title('Q1 - Section 4 - Transmit field picture');
    xlabel('X[mm]');ylabel('Z[mm]');
    
    % Convert Intensity results to dB
    a = (1-1e-4)/(max(Result(:))- min(Result(:)));
    b = 1 - a*max(Result(:));
    Result = a*Result + b;
    Result_dB=10*log10(Result);    
    subplot(1,2,2);
    imagesc(1000*x(:),1000*z(:),Result_dB);
    colormap(hot);
    title('Q1 - Section 4 - Transmit field picture dB');
    xlabel('X[mm]');ylabel('Z[mm]');
    colorbar;
    z_30mm = round(find(abs(z(:) - 0.03) == min(abs(z(:) - 0.03)),1)/100);
    z_60mm = round(find(abs(z(:) - 0.06) == min(abs(z(:) - 0.06)),1)/100);
    x = linspace(-10,10,100);
    figure('Name','Q1 - Section 4 - Literal intersects at z=30mm,60mm');
    plot(x, Result_dB(z_30mm,:), x,  Result_dB(z_60mm,:));
    xlabel('X[mm]')
    title('Q1 - Section 4 - Literal intersects at z=30mm,60mm');
    legend('z=30[mm]','z=60[mm]');
end
