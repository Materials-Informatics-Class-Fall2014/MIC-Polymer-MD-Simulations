% This is a general test code for trivial two-point statistics.
% The code first generates a trial microstructure given an assigned
% number of total local states and the number of elements per side 
% length (the code currently only works for 2-D microstructures).
% Then, two point statistics are performed for an autocorrelation
% and a crosscorrelation. The original microstructure, the 
% autocorrelation and the crosscorrelation are each plotted
% separately.

% Noah Paulson, September 16, 2014

close all

%% Microstructure Generation and Plotting

% H: number of local states
H = 2;

% el: side length of square
el = 21;

% S: the total number of spatial bins
S = el^2;

% select 0 for quilted microstructure or 1 for random microstructure
quilted_microstructure = 1;

if quilted_microstructure == 1
	% 'quilted' microstructure
	ms = zeros(el,el);
	for ii = 1:el
		for jj = 1:el
			inter = cos((2.5*ii)/(2*pi()))*cos((2.5*jj)/(2*pi()));
			ms(ii,jj) = round(1+(H-1)*(0.5 + 0.5*inter));
		end
	end
else
	% randomly generated microstructure    
	ms = round(1 + (H-1)*rand(el,el));
	% save('ms','ms')
	% mean(ms(:))
end

% plot the microstructure function
figure(1)
image(ms,'CDataMapping','scaled')
colormap('jet')
axis equal 
axis tight
colorbar

%% Calculate the microstructure function

% the real space microstructure function
msf = zeros(el,el,H);

for h = 1:H
    msf(:,:,h) = (ms == h);
end

% FFT frequency space microstructure function
MSF = fft2(msf);

%% 2-pt statistics autocorrelation

% 2-pt statistics autocorrelation for phase 1 in fft space
TPS_auto = (1 / (el^2)) .* abs(MSF(:,:,1)).^2;

% 2-pt statistics autocorrelation for phase 1 in real space
tps_auto = real(ifft2(TPS_auto));

% plot the autocorrelation

% center the autocorrelation for viewing
shi = floor(0.5*el);
% tps_auto_cen =  circshift(circshift(tps_auto,[shi,0]),[0,shi]); %Matlab 2013b version
tps_auto_cen =  circshift(circshift(tps_auto,shi,1),shi,2); %Matlab 2014b version

figure(2)
image(tps_auto_cen,'CDataMapping','scaled')
colormap('jet')
axis equal 
axis tight
title('2-Point statistics autocorrelation')
colorbar


%% 2-pt statistics crosscorrelation

% 2-pt statistics crosscorrelation for phase 1 in fft space
pha1 = 1;
pha2 = 2;
P1 = abs(MSF(:,:,pha1)) .* exp(-1i * angle(MSF(:,:,pha1)));
P2 = abs(MSF(:,:,pha2)) .* exp(1i * angle(MSF(:,:,pha2)));

TPS_cross = (1 / (el^2)) .* P1 .* P2;

% 2-pt statistics crosscorrelation for phase 1 in real space
tps_cross = real(ifft2(TPS_cross));

% plot the crosscorrelation

% center the crosscorrelation for viewing
shi = floor(0.5*el);
% tps_cross_cen =  circshift(circshift(tps_cross,[shi,0]),[0,shi]); %Matlab 2013b version
tps_cross_cen =  circshift(circshift(tps_cross,shi,1),shi,2); %Matlab 2014b version

figure(3)
image(tps_cross_cen,'CDataMapping','scaled')
colormap('jet')
axis equal 
axis tight
title('2-point statistics crosscorrelation')
colorbar
