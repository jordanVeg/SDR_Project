% SDR_Proj_4PAM_Tx.m
% Transmitter and channel impairment generator
% 4PAM encoding is used
% Output signal is "r"

%clear
% Enable plots of signals at various points as signal is developed
plotsONOFF=[0 0 0 0 0 0 0 0 0];
%plotsONOFF=[1 1 1 1 1 1 1 1 1];
% Enable transmitted signal impairments
noiseON      = 0;
multipathON  = 0;
freqErrorON  = 0;
phaseNoiseOn = 0;
BAUDerrorON  = 0;
%%

% Message Generation
%m='01234 I wish I were an Oscar Meyer wiener 56789';

m=['This is the first frame which you probably shouldn'  ...
   '''t be able to decode perfectly unless you "cheat" ' ...
   'and give your receiver the initial points.  Now we''re' ...
   ' into the fourth frame. You might be able to decod'  ...
   'e this one, and now the fifth.  The next frame is '  ...
   'where your receiver should definitely be operating'  ...
   'without errors.  You might want to re-test your re'  ...
   'ceiver by using different initial parameters, and '  ...
   'different stepsizes to see what the effect is.  It'  ...
   '''s probably helpful to plot the time history of th' ...
   'e adaptive parameter elements, too, so you can see'  ...
   ' if they''re taking too long to converge, if they s' ...
   'eem unstable, etc.  This easy test vector has no i'  ...
   'mpairments, i.e. no carrier frequency offset, a si'  ...
   'mple unity gain channel, no broadband noise, no na'  ...
   'rrowband interferers, no adjacent users, no baud t'  ...
   'iming clock offset, no phase noise.  So, the algor'  ...
   'ithms in your receiver should settle fairly quickl'  ...
   'y, and remain constant... If your receiver has mad'  ...
   'e it this far with no errors, and performs error-f'  ...
   'ree even when you change the initial parameter val'  ...
   'ues, Good job!!  Now it''s time to add impairments!   '];
%%

% Frame parameters
frameParams.userDataLength=100;  
frameParams.preamble='A0Oh well whatever Nevermind';     

% Channel parameters, Adjacent Users, Interference
if multipathON == 1
    chanParams.multipath=[0.5 1 -0.6];
else    
    chanParams.multipath=1;
end    
% chanParams.c1=[1 0 0];
% chanParams.c2=[1 0 0];
% chanParams.randomWalkVariance=0;
if noiseON == 1
    chanParams.SNR=30;  % dB
else
    chanParams.SNR=Inf;
end
% chanParams.adjacentUser1Power=-Inf;
% chanParams.adjacentUser1f_if=0;
% chanParams.adjacentUser1Chan=[1 0 0];
% chanParams.adjacentUser2Power=-Inf;
% chanParams.adjacentUser2f_if=0;
% chanParams.adjacentUser2Chan=[1 0 0];
% chanParams.NBIfreq=0;
% chanParams.NBIPower=-Inf;

% RF Parameters
if freqErrorON==1
    txParams.f_if_err=1;  % Percentage Error in IF frequency
else
    txParams.f_if_err=0;  % No IF frequency error
end    
if BAUDerrorON==1
    txParams.T_t_err=0.1;   % Percentage Error in Baud rate (symbol rate)
else
    txParams.T_t_err=0;   % No Baud rate error
end    
if phaseNoiseOn == 1
    txParams.phaseNoiseVariance=0.0001;  % Add phase noise
else
    txParams.phaseNoiseVariance=0;  % No phase noise
end
rfParams.SRRCLength=4;  % Square-Root Raised Cosine filter length
rfParams.SRRCrolloff=0.3; % SRRC rolloff factor, Beta
rfParams.f_s=8400;       % Sampling frequency
rfParams.T_t=0.0064;    % nominal symbol period
rfParams.f_if=2000;     % intermediate frequency of the sampled input 
                        % to the digital receiver

%%
% TRANSMITTER
% function [r, s]=BigTransmitter(m, frameParams, rfParams, chanParams)

% DIVIDE DATA INTO FRAMES
% re-organize m so that each row is a frame
linesOfText=floor(size(m,2)/frameParams.userDataLength);
m2=reshape(m(1,1:frameParams.userDataLength*linesOfText),...
    frameParams.userDataLength,linesOfText)';
%%
% ADD FRAME MARKERS
% insert frame marker (preamble)
m3=reshape([repmat([frameParams.preamble],linesOfText,1) m2(:,:)]',...
    (length(frameParams.preamble)+...
    frameParams.userDataLength)*linesOfText,1);
%%
% ENCODE MESSAGE TO SYMBOLS
% encode characters into 4-PAM source vector
s=letters2pam(m3);

%%

% generate received signal 
% (transmitter RF frontend --> channel --> sampled-IF receiver)
% r0=Tx_rf(s, rfParams, chanParams);

%%
% CONVERT SYMBOLS INTO ANALOG PULSES

% introduce error into transmitter IF frequency and baud timing
txParams.f_if_tx=rfParams.f_if*(1+txParams.f_if_err/100);
txParams.T_t_tx=rfParams.T_t*(1+txParams.T_t_err/100);

% determine suitable oversampling/downsampling factor
rate_ratio=rfParams.f_s*txParams.T_t_tx;
[M,N]=rat(rate_ratio);
    % We want to have the output of this transmitter, channel, and receiver analog
    % front end to be sampled at the f_s rate. But the T_t_tx symbol period
    % is such that the ratio of the sampling rate to symbol rate is not an
    % integer value.  So, we oversample by M, do the pulse-shaping, carrier
    % modulation, etc. and then downsample by N to get to the desired
    % sampling rate.
    % e.g. if f_s=850 and T_t_tx=0.0064, the rate ratio is f_s*T_t_tx=5.44
    %      so we want M/N=5.44 where M and N are integers.
    %      With M=136 and N=25, M/N=5.44

disp(['Sampling Rate = ', num2str(rfParams.f_s)])
disp(['Symbol Rate   = ', num2str(1/txParams.T_t_tx)])
disp(['ratio of Sampling Rate to Symbol Rate: ', num2str(rate_ratio)])    
disp('oversampling/downsampling factor: ');
disp(['     M = ', num2str(M),', N = ', num2str(N)]);
disp(['Oversampling Rate = ', num2str(M*1/txParams.T_t_tx)])
disp(['Carrier Frequency = ', num2str(rfParams.f_if)])

% generate pulse-shaped signal
sm=upsample(s,M);    % upsample by M
psf=srrc(rfParams.SRRCLength,rfParams.SRRCrolloff,M,0); % pulse-shaping filter
x=conv(psf,sm);      % pulse shaping of symbols

% create time vector, t
Ts=txParams.T_t_tx/M;
t=0:Ts:Ts*(length(x)-1);

tRange=floor(length(t)/2):floor(length(t)/2)+25*M;  % t range for plot

fig=1; % Plot encoded message
if plotsONOFF(fig)==1
    figure(fig)
    square_s=repelem(s,M);     % 4-level signal for plotting
    plotspec4(square_s,Ts,t,tRange), zoom xon
    subplot(2,1,1),title('Message as 4-PAM symbols')
    clearvars square_s
end

fig=1; % Plot pulse-shaping filter response
if plotsONOFF(fig)==1
    figure(fig)
    freqz(psf,512,0:(pi/8)/512:pi/8), zoom xon
    subplot(2,1,1),title('Pulse-Shaping filter response')
end

tRange=(floor(length(t)/2):floor(length(t)/2)+25*M)+4*M;  % t range for plot

fig=fig+1; % Plot pulse-shaping filter response
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(x,Ts,t,tRange), zoom xon, axis([-1/(Ts*M) 1/(Ts*M),ylim])
    subplot(2,1,1),title('Pulse-Shaped Message Signal')
end

clear sm
%%

% ADD PHASE NOISE TO CARRIER
if phaseNoiseOn == 1
    % generate phase noise process
    p_noise=cumsum(randn(size(x))*sqrt(txParams.phaseNoiseVariance/N));
    % create carrier with phase noise
    carrier=cos(2*pi*txParams.f_if_tx*t+p_noise); 
    
    fig=fig+1;  % plot phase noise
    if plotsONOFF(fig)==1
        figure(fig)
        plotspec4(p_noise,Ts,t,tRange), zoom xon
        subplot(2,1,1),title('Carrier Phase Noise Signal')
    end
else 
    carrier=cos(2*pi*txParams.f_if_tx*t); % no phase noise
end

% plot carrier
tRange=floor(length(t)/2):floor(length(t)/2)+M;
fLim=min([2*rfParams.f_if,0.5*M/txParams.T_t_tx]);
fig=fig+1;
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(carrier,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1),title('Transmitter Carrier Signal')
end

clearvars p_noise      % free up a little memory 

% modulate carrier with signal to produce RF signal (well, technically IF)
x_rf=carrier.*x;

clear x carrier        % free up a little memory 

tRange=floor(length(t)/2):floor(length(t)/2)+25*M;


fig=fig+1;  % plot transmitted 4-PAM signal
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(x_rf,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1),title('Transmitted 4-PAM Signal')
end
%%

% ADD CHANNEL IMPAIRMENTS
% add multipath to RF signal
r=filter(upsample(chanParams.multipath,M),1,x_rf);  
                % needed to upsample multipath before convolving with signal

clearvars x_rf

fig=fig+1;  % plot signal plus multipath
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(r,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1),title('Signal with Multipath')
end

% add noise to RF signal
sp=pow(r);                      
np=sp/10^(chanParams.SNR/10);
n=sqrt(np)*randn(1,length(r));
r=r+n;

clearvars n

fig=fig+1;  % plot signal plus noise
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(r,Ts,t,tRange), zoom xon, axis([1800,2200,ylim])
    subplot(2,1,1),title('Signal Plus Noise')
end

% perform downsampling, normalization
% this is the output of the sampler at the receiver
r=r(1:N:end);
r=r*sqrt(length(r))/norm(r); % r = OUTPUT OF SAMPLER AT RECEIVER

Ts=Ts*N;
t=0:Ts:Ts*(length(r)-1);
tRange=floor(length(t)/2):floor(length(t)/2)+25*floor(M/N);


fig=fig+1;  % plot output of sampler at receiver
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(r,Ts,t,tRange), zoom xon, axis([1800,2200,ylim])
    subplot(2,1,1),title('Output of Receiver Sampler')
end