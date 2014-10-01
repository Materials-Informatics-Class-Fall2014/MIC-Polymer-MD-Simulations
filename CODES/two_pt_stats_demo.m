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
el = 11;

% S: the total number of spatial bins
S = el^2;

% select 0 for quilted microstructure or 1 for random microstructure
quilted_microstructure = 1;

if quilted_microstructure == 1
	% 'quilted' microstructure
	ms = zeros(el,el);
	for ii = 1:el
		for jj = 1:el
			inter = cos((1.25*ii)/(2*pi()))*cos((1.25*jj)/(2*pi()));
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


%% Compute the autocorrelations

% find the volume fraction of phase 1
msf1 = msf(:,:,1);
% vol1: volume fraction of phase 1
vol1 = sum(msf1(:))/S;

TPS = zeros(el,el,H,H);
tps = zeros(el,el,H,H);

for h = 1: H
    
    % compute the correlation
    TPS(:,:,h,h) = (1/S)*(abs(MSF(:,:,h)).^2);
    % convert the correlation to real space
    tps(:,:,h,h) = real(ifft2(TPS(:,:,h,h)));

end


%% 2-pt statistics autocorrelation and crosscorrelations (2 through H-1) for phase 1

for h=2: H - 1
    
    P1 = abs(MSF(:,:,1)) .* exp(-1i * angle(MSF(:,:,1)));
    P2 = abs(MSF(:,:,h)) .* exp(1i * angle(MSF(:,:,h)));

    % compute the correlation
    TPS(:,:,1,h) = (1 / (el^2)) .* P1 .* P2;

    % correlation in real space
    tps(:,:,1,h) = real(ifft2(TPS(:,:,1,h)));

end

%% compute the F1H crosscorrelation

tps(:,:,1,H) = vol1 .* ones(el,el);
for h = 1:H-1
    tps(:,:,1,H) = tps(:,:,1,H) - tps(:,:,1,h);
end

TPS(:,:,1,H) = fft2(tps(:,:,1,H));


%% compute the rest of the statistics

for ii = 2:H
    for jj = 1:H
        
        if ii ~= jj
            TPS(:,:,ii,jj) = (conj(TPS(:,:,1,ii)).*TPS(:,:,1,jj)) ./ TPS(:,:,1,1);
            tps(:,:,ii,jj) = real(ifft2(TPS(:,:,ii,jj)));
        end
        
    end
end


%% plot all of the 2-pt statistics

shi = floor(0.5*el);
c = 1;
tps_cen = zeros(el,el,H,H);

figure(2)

for ii = 1:H
    for jj = 1:H
    
    % center the crosscorrelation for viewing

    tps_cen(:,:,ii,jj) =  circshift(circshift(tps(:,:,ii,jj),[shi,0]),[0,shi]); %Matlab 2013b version
%     tps_cen(:,:,ii,jj) =  circshift(circshift(tps(:,:,ii,jj),shi,1),shi,2); %Matlab 2014b version
               
    subplot(H,H,c)
    image(tps_cen(:,:,ii,jj),'CDataMapping','scaled')
    colormap('jet')
    axis equal 
    axis tight
    title(['2-point statistics: ',num2str(ii),' ,',num2str(jj)])
    colorbar
    
    c = c + 1;
    end
end


%% Test the 2-pt statistics: (less computationaly efficient)
% 
% figure(3)
% 
% TPS_t = zeros(el,el,H,H);
% tps_t = zeros(el,el,H,H);
% tps_t_cen = zeros(el,el,H,H);
% 
% c = 1;
% 
% for ii= 1:H
%     for jj = 1:H
%     
%     P1 = abs(MSF(:,:,ii)) .* exp(-1i * angle(MSF(:,:,ii)));
%     P2 = abs(MSF(:,:,jj)) .* exp(1i * angle(MSF(:,:,jj)));
% 
%     % compute the correlation
%     TPS_t(:,:,ii,jj) = (1 / (el^2)) .* P1 .* P2;
% 
%     % correlation in real space
%     tps_t(:,:,ii,jj) = real(ifft2(TPS_t(:,:,ii,jj)));
% 
%     tps_t_cen(:,:,ii,jj) =  circshift(circshift(tps_t(:,:,ii,jj),shi,1),shi,2); %Matlab 2014b version
%                
%     subplot(H,H,c)
%     image(tps_cen(:,:,ii,jj),'CDataMapping','scaled')
%     colormap('jet')
%     axis equal 
%     axis tight
%     title(['2-point statistics test: ',num2str(ii),' ,',num2str(jj)])
%     colorbar
%     
%     c = c + 1;
%     
%     end
% end
