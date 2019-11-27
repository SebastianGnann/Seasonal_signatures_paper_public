function plotFourierSpectrum(F,Q,w,ID,title_str)
%PLOTFOURIERSPECTRUM 
%
% INPUT
% F: forcing, typically approximated by P - PET [mm]
% Q: streamflow [mm]
% w: angular frequency [1/days]
% ID: catchment ID
% title: catchment title
%
% OUTPUT 
% plot
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk


%% Fourier transform
N = length(Q(:,2)); % length of time series
fs = 1; % 1d sample frequency
% f = (0:N-1)*fs/N; % frequency range 1/d

% forcing
fft_F = fft(F(:,2)); % compute Fourier transform
P2_F = abs(fft_F/N); % compute Fourier transform
P1_F = P2_F(2:floor(N/2)+1); % compute one-sided spectrum
P1_F(2:end-1) = 2*P1_F(2:end-1);

% streamflow
fft_Q = fft(Q(:,2)); % compute Fourier transform
P2_Q = abs(fft_Q/N); % compute Fourier transform
P1_Q = P2_Q(2:floor(N/2)+1); % compute one-sided spectrum
P1_Q(2:end-1) = 2*P1_Q(2:end-1);

f = fs*(1:floor(N/2))/N; % define frequency domain
t = 1./f/365; % and asociated time domain in years

% plot
n = 10;
colour_mat = (brewermap(n,'Set1')); %flip
fig1 = figure('Name',title_str,'NumberTitle','off','pos',[10 10 300 250]); % 

subplot(2,1,1)
title(title_str)
hold on
plot(t,P1_F,'Linewidth',1.5,'color',colour_mat(2,:));
% xlabel('Time [years]')
ylabel('|P1(P-E_p)|')
set(gca,'xscale','log')
xtick=[1/100 1/10 1.0 10. 100.];
xticklab = {'0.01','0.1','1','10','100'};
set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex')

subplot(2,1,2)
hold on
plot(t,P1_Q,'Linewidth',1.5,'color',colour_mat(5,:));
xlabel('Time [years]')
ylabel('|P1(Q)|')
set(gca,'xscale','log')
xtick=[1/100 1/10 1.0 10. 100.];
xticklab = {'0.01','0.1','1','10','100'};
set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex')

% save fig
set(fig1,'Units','Inches');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3),pos(4)]);
fig_name_raw = strcat('Fourier_spectrum_','_',num2str(ID));
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(fig1,fig_path,'-dpdf','-r500');

end

