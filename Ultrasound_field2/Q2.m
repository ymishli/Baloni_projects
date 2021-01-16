clc; close all;clear all;

% Initialize Field
%addpath('C:\Users\Rony\Desktop\Ultrasound\Field2\Field_II_ver_3_22_windows');
field_init(0)

% Control the plots plot_on = [sec_a sec_b sec_c sec_d sec_e]
plot_on = [0 1 1 1 1];

% Generate the transducer aperture for send and receive
f0          = 2.5e6;            %  Transducer center frequency [Hz]
% f0          = 5e6;            %  Transducer center frequency [Hz]
fs          = 100e6;            %  Sampling frequency [Hz]
c           = 1490;             %  Speed of sound [m/s]
lambda      = c/f0;             %  Wavelength [m]
width       = 0.29/1000;        %  Width of element
height      = 1/1000;          % Height of element [m]
kerf        = 0.010/1000;                % Kerf [m]
focus       = [0 0 60]/1000;    % Fixed focal point [m]
N_elements  = 64;                % Number of physical elements
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
    figure('Name','Q2 - Section 1 - Transducer Geometry');
    show_xdc_geir(tx,1);    
    axis equal;
    view(3);
    zlim([-5 5]);
    title('Q2 - Section 1 - Transducer Geometry');
end

%1.b + 1.c
Bw = 0.6;
t_h = (-2/f0:1/fs:2/f0);
impulse_response = gauspuls(t_h,f0,Bw);    
xdc_impulse(tx,impulse_response);    
for field_num=1:3 % I run it three times for all the cases in section 2 and 3
    if field_num == 1
        xdc_focus(tx, 0, [0 0 40]/1000);
        original_focus_data = xdc_get(tx,'focus');        
        time_taps_generated_by_xdc_focus = original_focus_data(2:end);
    else
        pitch = width + kerf;
        plot_for_debug = 1;
        if field_num == 2
            N_elements_delay = Delay(N_elements,c,pitch,[0 0 40]/1000, plot_for_debug, time_taps_generated_by_xdc_focus);
        end
        if field_num == 3
            N_elements_delay = Delay(N_elements,c,pitch,[5 0 40]/1000, plot_for_debug, time_taps_generated_by_xdc_focus);
        end
        N_elements_delay_matrix = repmat(N_elements_delay,N_elements,1);
        xdc_focus_times(tx,zeros(length(N_elements_delay),1),N_elements_delay_matrix);
    end
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
    if plot_on(2) == 1 || plot_on(3) == 1
        if field_num == 1
            figure('Name','Q2 - Section 2 - Transmit field picture at [-10,10]mm*[10,80]mm');
        elseif field_num == 2
            figure('Name','Q2 - Section 3 - Transmit field with our Delay taps focus at [0 0 40]');
        else % field_num == 3
            figure('Name','Q2 - Section 3 - Transmit field with our Delay taps focus at [5 0 40]');
        end                
        subplot(1,2,1);
        imagesc(1000*x(:),1000*z(:),Result_norm);
        colormap(hot)
        if field_num == 1
            title('Q2 - Section 2 - Transmit field picture');
        elseif field_num == 2
            title('Q2 - Section 3 - Transmit field with our Delay taps focus at [0 0 40]');
        else % field_num == 3
            title('Q2 - Section 3 - Transmit field with our Delay taps focus at [5 0 40]');
        end                        
        xlabel('X[mm]');ylabel('Z[mm]');
        % Convert Intensity results to dB
        a = (1-1e-4)/(max(Result(:))- min(Result(:)));
        b = 1 - a*max(Result(:));
        Result = a*Result + b;
        Result_dB=10*log10(Result);            
        subplot(1,2,2);
        imagesc(1000*x(:),1000*z(:),Result_dB);
        colormap(hot);
        title('Q2 - Section 2 - Transmit field picture dB');
        xlabel('X[mm]');ylabel('Z[mm]');
        colorbar;
        if field_num == 1            
            x = linspace(-10,10,100);
            z_20mm = round(find(abs(z(:) - 0.02) == min(abs(z(:) - 0.02)),1)/100);
            z_40mm = round(find(abs(z(:) - 0.04) == min(abs(z(:) - 0.04)),1)/100);
            z_40_cross_section = Result_dB(z_40mm,:);
            figure('Name','Q2 - Section 2 - Literal cross sections at z=20mm,40mm');
            plot(x, Result_dB(z_20mm,:), x,  Result_dB(z_40mm,:));            
            xlabel('X[mm]')
            title('Q2 - Section 2 - Literal cross sections at z=20mm,40mm');
            legend('z=20[mm]','z=40[mm]');
        end        
    end
end

%1.d
xdc_focus(tx, 0, [0 0 40]/1000);
[x,y,z]=meshgrid(linspace(-10,10,100)/1000,0,linspace(10,80,100)/1000);
points=[x(:) y(:) z(:)];
for type_of_window = 1:3 % 1 - rect , 2 - guasian , 3 - cosian :-)
    if type_of_window == 1 % rect
        xdc_apodization(tx,zeros(N_elements,1),ones(N_elements,N_elements));
    elseif type_of_window == 2 % guasian
        guasian_apo = hamming(N_elements);
        guasian_apo_matrix = repmat(guasian_apo' ,N_elements,1);
        xdc_apodization(tx,zeros(N_elements,1),guasian_apo_matrix);
    else % cosian
        cos_apo = tukeywin(N_elements,1);
        cos_apo_matrix = repmat(cos_apo' ,N_elements,1);
        xdc_apodization(tx,zeros(N_elements,1),cos_apo_matrix);
    end    
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
    a = (1-1e-4)/(max(Result(:))- min(Result(:)));
    b = 1 - a*max(Result(:));
    Result = a*Result + b;
    Result_dB=10*log10(Result);
    z_40mm = round(find(abs(z(:) - 0.04) == min(abs(z(:) - 0.04)),1)/100);    
    if plot_on(4) == 1
        if type_of_window == 1 % open figure in first itertion only
            figure('Name','Q2 - Section 4 - Apodization effect on z=40[mm] cross section');            
        end    
        plot(linspace(-10,10,100),  Result_dB(z_40mm,:));
        hold on;
        if type_of_window == 3 % do cosmetic stuff on graph for last iteration
            xlabel('X[mm]');ylabel('Pressure');
            title('Q2 - Section 4 - Apodization effect on z=40[mm] cross section');
            legend('rect','guasian','cos');
            hold off;
        end    
    end    
end
   
%1.e
% return the apodization to rect()
xdc_apodization(tx,zeros(N_elements,1),ones(N_elements,N_elements));
% Generate aperture for receive
rx = xdc_linear_array (N_elements, width, height, kerf, N_sub_x, N_sub_y, focus);
% Set excitation
t = (0:1/fs:1.5/f0);
excitaion = sin(2*pi*f0*t);
xdc_excitation(rx,excitaion());
% Set impulse response
Bw = 0.6;
t_h = (-2/f0:1/fs:2/f0);
impulse_response = gauspuls(t_h,f0,Bw);    
xdc_impulse(rx,impulse_response);    
% Define points to calc the received signal from
[x,y,z]=meshgrid(linspace(-10,10,100)/1000,0,linspace(10,80,100)/1000);
points=[x(:) y(:) z(:)];
Im_size=[length(x),1,length(z)];
[hp,start_t]=calc_hhp(tx,rx,points);
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
a = (1-1e-4)/(max(Result(:))- min(Result(:)));
b = 1 - a*max(Result(:));
Result = a*Result + b;
Result_dB=10*log10(Result);
z_40mm = round(find(abs(z(:) - 0.04) == min(abs(z(:) - 0.04)),1)/100);
if plot_on(5) == 1
    figure('Name','Q2 - Section 5 - Received vs Transmit z=40[mm] cross section');
    plot(linspace(-10,10,100), Result_dB(z_40mm,:), linspace(-10,10,100), z_40_cross_section);
    xlabel('X[mm]'); ylabel('Pressure');
    title('Q2 - Section 5 - Received vs Transmit z=40[mm] cross section');
    legend('Received','Transmitted');
end

%1.f
grating_lobes_45_deg = sqrt(2)*N_elements*c/f0 - width;
% kerf_vec        = (0.010:0.010:0.090)/1000;                % Kerf [m]
% kerf_vec        = linspace(kerf,grating_lobes_45_deg, 6);
kerf_vec        = linspace(kerf,80*kerf, 6);
[x,y,z]=meshgrid(linspace(-20,20,100)/1000,0,linspace(10,80,100)/1000);
points=[x(:) y(:) z(:)];
Im_size=[length(x),1,length(z)];
% Generate aperture for transmission
figure('Name','Q2 - Section 6 - Detecing Grating Lobes with kerf increase');
for kerf_index = 1:length(kerf_vec)
    tx = xdc_linear_array (N_elements, width, height, kerf_vec(kerf_index), N_sub_x, N_sub_y, focus);

    % Set the excitation of the transmit aperture
    t = (0:1/fs:1.5/f0);
    excitaion = sin(2*pi*f0*t);
    xdc_excitation(tx,excitaion());
    %1.b + 1.c
    Bw = 0.6;
    t_h = (-2/f0:1/fs:2/f0);
    impulse_response = gauspuls(t_h,f0,Bw);    
    xdc_impulse(tx,impulse_response);
    xdc_focus(tx, 0, [0 0 40]/1000);
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
    a = (1-1e-4)/(max(Result(:))- min(Result(:)));
    b = 1 - a*max(Result(:));
    Result = a*Result + b;
    Result_dB=10*log10(Result);
    subplot(2,3,kerf_index);
    imagesc(1000*x(:),1000*z(:),Result_norm);    
    xlabel('X[mm]');ylabel('Z[mm]');
    colormap(hot)
    subtitle = strcat('kerf = ', num2str(kerf_vec(kerf_index)));
    title(subtitle);
end
