close all
clear all
clc

% Initialize Field
field_init(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field code for simulating basic linear arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0=3e6;                     %  Transducer center frequency [Hz]
fs=100e6;                   %  Sampling frequency [Hz]
c=1540;                     %  Speed of sound [m/s]
lambda=c/f0;                %  Wavelength [m]
width=lambda;               %  Width of element
element_height=5/1000;      %  Height of element [m]
kerf=0.1/1000;              %  Kerf [m]
focus=[0 0 70]/1000;        %  Fixed focal point [m]
N_tx_elements=128;          %  Number of physical elements in the transmit aperture
N_rx_elements=128;          %  Number of physical elements in the receive aperture

%  Set the relevent simulation parameters
set_sampling(fs);                   %  Sets sampling frequency
set_field('use_triangles',0);       %  Tells whether to use triangles (1) or not (0)
set_field('use_rectangles',1);      %  Tells whether to use rectangles (1) or not (0)
set_field('use_att',0);             %  Tells whether to use attenuation (1) or not (0)
set_field('c',c);                   %  Sets the speed of sound

%  Generate aperture for transmission
tx=xdc_linear_array(N_tx_elements,width,element_height,kerf,1,1,focus);

%  Set the impulse response of the transmit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

xdc_impulse(tx,impulse_response);

%  Set the excitation of the transmit aperture
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(tx,excitation);

figure;
subplot(2,1,1);
plot(0:1/fs:2/f0,excitation);title('Excitation Signal');xlabel('t[sec]');
subplot(2,1,2);
plot(0:1/fs:2/f0,impulse_response);title('Impulse Response');xlabel('t[sec]');


%  Generate aperture for reception
rx=xdc_linear_array(N_rx_elements,width,element_height,kerf,1,1,focus);

%  Set the impulse response for the receive aperture
xdc_impulse(rx,impulse_response);

% Apodize the array on transmit and/or receive
xdc_apodization(tx,0,ones(1,N_tx_elements));
xdc_apodization(rx,0,ones(1,N_rx_elements));

% xdc_apodization(tx,0,hanning(N_tx_elements)');
% xdc_apodization(rx,0,hanning(N_rx_elements)');

% Define matrix of [x y z] coordinates for where the field calculations occur
W=5/1000;%[m]
[x,y,z]=meshgrid(-W/2:0.0001:W/2,0,0);
z=z+focus(3);
points=[x(:) y(:) z(:)];

% Calculate the spatial impulse response for the transmit aperture
[h_tx,start_time_h_tx]=calc_h(tx,points);

% Calculate the spatial impulse response for the receive aperture (if different than tx)
[h_rx,start_time_h_rx]=calc_h(rx,points);

% Calculate the emitted field
[hp,start_time_hp]=calc_hp(tx,points);

% Calculate the pulse echo field
[hhp,start_time_hhp]=calc_hhp(tx,rx,points);

figure;

% Display the spatial impulse response for the transmit aperture
lateral_dim=1000*x';%[mm]
axial_dim=1000*c*(start_time_h_tx+1/fs*(0:size(h_tx,1)-1));
subplot(2,3,1);imagesc(lateral_dim,axial_dim,h_tx);colormap(gray);
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
title('Tx spatial impulse response');

% Display the spatial impulse response for the receive aperture
lateral_dim=1000*x';%[mm]
axial_dim=1000*c*(start_time_h_rx+1/fs*(0:size(h_rx,1)-1));
subplot(2,3,2);imagesc(lateral_dim,axial_dim,h_rx);colormap(gray);
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
title('Rx spatial impulse response');

% Display the emitted field
lateral_dim=1000*x';%[mm]
axial_dim=1000*c*(start_time_hp+1/fs*(0:size(hp,1)-1));
subplot(2,3,3);imagesc(lateral_dim,axial_dim,hp);colormap(gray);
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
title('Transmitted field');

% Display the pulse echo field (also PSF)
lateral_dim=1000*x';
axial_dim=1000*c/2*(start_time_hhp+1/fs*(0:size(hhp,1)-1));
subplot(2,3,4);imagesc(lateral_dim,axial_dim,hhp);colormap(gray);
xlabel('Lateral distance (mm)');
ylabel('Depth (mm)');
title('Pulse echo field (PSF)');

% Envelope detect the PSF
hhp_det=abs(hilbert(hhp));

lateral_dim=1000*x';
axial_dim=1000*c/2*(start_time_hhp+1/fs*(0:size(hhp,1)-1));
center_line=(size(hhp,2)-1)/2+1;

% Normalize the PSF and plot the axial cross-section
hhp_norm=hhp/max(max(hhp_det));
subplot(2,3,5);plot(axial_dim,hhp_norm(:,center_line));axis tight;
xlabel('Depth (mm)');
ylabel('Normalized Pressure');
title('Axial cross-section');

% Normalize the detected PSF and plot the lateral cross-section
hhp_det_norm_dB=20*log10(hhp_det/max(max(hhp_det)));
ind=rem(find(hhp_det==max(max(hhp_det))),size(hhp,1));
subplot(2,3,6);plot(lateral_dim,hhp_det_norm_dB(ind,:));axis tight;
xlabel('Lateral distance (mm)');
ylabel('Normalized Pressure (dB)');
title('Lateral cross-section');

%% Display transmission

W=10/1000;%[m]
[x,y,z]=meshgrid(-W/2:0.0001:W/2,0,-W:0.0001:W);
Im_size=[length(-W/2:0.0001:W/2),1,length(-W:0.0001:W)];
z=z+focus(3);
points=[x(:) y(:) z(:)];
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
 
% Close Field
field_end;

 