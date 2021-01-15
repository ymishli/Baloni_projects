close all
clear all
clc

% Initialize Field
field_init(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field code for simulating basic linear arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0=2.5e6;                     %  Transducer center frequency [Hz]
fs=100e6;                   %  Sampling frequency [Hz]
c=1490;                     %  Speed of sound [m/s]
lambda=c/f0;                %  Wavelength [m]
width       = 18.5/1000;        %  Width of element
element_height      = 13/1000;          % Height of element [m]
kerf=0;              %  Kerf [m]
focus=[0 0 60]/1000;        %  Fixed focal point [m]
N_tx_elements=1;          %  Number of physical elements in the transmit aperture

%  Set the relevent simulation parameters
set_sampling(fs);                   %  Sets sampling frequency
set_field('use_triangles',0);       %  Tells whether to use triangles (1) or not (0)
set_field('use_rectangles',1);      %  Tells whether to use rectangles (1) or not (0)
set_field('use_att',0);             %  Tells whether to use attenuation (1) or not (0)
set_field('c',c);                   %  Sets the speed of sound

%  Generate aperture for transmission
tx=xdc_linear_array(N_tx_elements,width,element_height,kerf,1,1,focus);

%  Set the impulse response of the transmit aperture
Bw = 0.6;
t_h = (-2/f0:1/fs:2/f0);
impulse_response = gauspuls(t_h,f0,Bw);
impulse_response = impulse_response.* sin(2*pi*f0*t_h);
xdc_impulse(tx,impulse_response);


%  Set the excitation of the transmit aperture
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(tx,excitation);

%% Display transmission

W=20/1000;%[m]
Z1 = -50/1000;
Z2 = 20/1000;
% x = linspace(-10,10,100);
% z = linspace(10,80,100);
% [x,y,z]=meshgrid(-W/2:0.0001:W/2,0,-W:0.0001:W);
% [x,y,z]=meshgrid(-W/2:0.0002:W/2,0,-W:0.0002:W);
[x,y,z]=meshgrid(-W/2:0.0002:W/2,0,Z1:0.0007:Z2);
% [x,y,z]=meshgrid(x,0,z);
% Im_size=[length(x),1,length(z)];
% Im_size=[length(-W/2:0.0002:W/2),1,length(-W:0.0002:W)];
Im_size=[length(-W/2:0.0002:W/2),1,length(Z1:0.0007:Z2)];
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
 a = (1-1e-4)/(max(Result(:))- min(Result(:)));
 b = 1 - a*max(Result(:));
 Result = a*Result + b;
 Result_dB=10*log10(Result);
 subplot(1,2,2);
 imagesc(1000*x(:),1000*z(:),Result_dB);
 colormap(hot);
 title('Transmit Field Example dB');
 xlabel('X[mm]');ylabel('Z[mm]');
 colorbar
 
% Close Field
field_end;

 