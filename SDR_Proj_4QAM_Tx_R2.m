% SDR_Proj_4QAM_Tx.m
% Transmitter and channel impairment generator
% 4QAM encoding is used
% Output signal is "r"
%% 
% Enable plots of signals at various points as the received signal is developed
% 
% Set bits to 0 = plot off or 1 = plot on

plotsONOFF=[0];  % pulse-shaping filter frequency response
plotsONOFF=[plotsONOFF,0];  % 1/2 of Message as 2-PAM Symbols
plotsONOFF=[plotsONOFF,0];  % other 1/2 of Message as 2-PAM Symbols
plotsONOFF=[plotsONOFF,0];  % Pulse-Shaped 1/2 Message Signal
plotsONOFF=[plotsONOFF,0];  % Phase Noise Signal Added to Carriers
plotsONOFF=[plotsONOFF,0];  % Transmitter Cosine Carrier Signal
plotsONOFF=[plotsONOFF,0];  % Transmitter QAM Signal
plotsONOFF=[plotsONOFF,0];  % Multipath Impulse Response
plotsONOFF=[plotsONOFF,0];  % Signal with Multipath
plotsONOFF=[plotsONOFF,0];  % Signal Plus Noise
plotsONOFF=[plotsONOFF,1];  % Received and Downsampled QAM Signal

%% 
% Enable transmitted signal impairments

noiseON=0;
multipathON = 0;
freqErrorON = 0;
phaseNoiseON = 0;
BAUDerrorON = 0;

phaseErrorOn = 0; % not using this one (fixed phase offset)
%% 
% Message Generation

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
% Define parameters that characterize transmitted signal Frame parameters

frameParams.userDataLength=100; %length(m); 
frameParams.preamble='A0Oh well whatever Nevermind';     

% Channel parameters, Adjacent Users, Interference
if multipathON == 1
    chanParams.multipath=[0.1 1 -0.1];
else    
    chanParams.multipath=1;
end    
% chanParams.c1=[1 0 0];
% chanParams.c2=[1 0 0];
% chanParams.randomWalkVariance=0;
if noiseON == 1
    chanParams.SNR=30;
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
    txParams.f_if_err=0;  % Percentage Error in IF frequency
end    
if BAUDerrorON==1
    txParams.T_t_err=1;   % Percentage Error in Baud rate (symbol rate)
else
    txParams.T_t_err=0;   % Percentage Error in Baud rate (symbol rate)
end    
if phaseNoiseON == 1
    txParams.phaseNoiseVariance=0.0001;  % No phase noise
else
    txParams.phaseNoiseVariance=0;       % No phase noise
end
rfParams.SRRCLength=4;    % Square-Root Raised Cosine filter length
rfParams.SRRCrolloff=0.3; % SRRC rolloff factor, Beta
rfParams.f_s=8400;        % Sampling frequency
rfParams.T_t=0.0064;      % Symbol period
rfParams.f_if=2000;       % IF frequency

disp('PARAMTERS USED IN THIS SIMULATION'),disp(' ')
disp('Tx Parameters not known to the Rx'),disp(txParams)
disp('Channel Parameters not known to the Rx'),disp(chanParams)
disp('Frame Parameters known to both Tx and Rx'),disp(frameParams)
disp('RF Parameters known to both the Tx and Rx'),disp(rfParams)

% function [r, s]=BigTransmitter(m, frameParams, rfParams, chanParams)
% Brought in the code (below) needed to implement the parts of BigTransmitter
% needed.
%% 
% Divide message into two messages using every other character

ms=m(1:2:end); % message for sine carrier
mc=m(2:2:end); % message for cosine carrier
%% 
% determine number of frames needed

numFRAMESms=ceil(size(ms,2)/frameParams.userDataLength);
numFRAMESmc=ceil(size(mc,2)/frameParams.userDataLength);
%% 
% pad messages with space characters

ms=[ms,repmat(' ',[1,numFRAMESms*frameParams.userDataLength-length(ms)])];
mc=[mc,repmat(' ',[1,numFRAMESms*frameParams.userDataLength-length(mc)])];
%% 
% reshape matrix so that rows are frames

ms2=reshape(ms(1,1:frameParams.userDataLength*numFRAMESms),...
    frameParams.userDataLength,numFRAMESms)';
mc2=reshape(mc(1,1:frameParams.userDataLength*numFRAMESmc),...
    frameParams.userDataLength,numFRAMESmc)';
%% 
% insert preamble at beginning of each frame

ms3=reshape([repmat([frameParams.preamble],numFRAMESms,1) ms2(:,:)]',...
    (length(frameParams.preamble)+...
    frameParams.userDataLength)*numFRAMESms,1);
mc3=reshape([repmat([frameParams.preamble],numFRAMESmc,1) mc2(:,:)]',...
    (length(frameParams.preamble)+...
    frameParams.userDataLength)*numFRAMESmc,1);
%% 
% encode characters into binary source vectors

ms4=zeros(1,8*length(ms3));            % storage for ascii 8-bit binary code
for k=0:length(ms3)-1                  % code characters in 8-bit ascii
  bc=dec2base(double(ms3(k+1)),2,8);
  for i=0:7
    ms4(8*k+1+i)=2*bin2dec(bc(i+1))-1;
  end  
end
mc4=zeros(1,8*length(mc3));            % storage for ascii 8-bit binary code
for k=0:length(mc3)-1                  % code characters in 8-bit ascii
  bc=dec2base(double(mc3(k+1)),2,8);
  for i=0:7
    mc4(8*k+1+i)=2*bin2dec(bc(i+1))-1;
  end  
end
%% 
% generate received signal (transmitter RF frontend --> channel --> sampled-IF 
% receiver)

% introduce error into transmitter IF and baud timing
txParams.f_if_tx=rfParams.f_if*(1+txParams.f_if_err/100);
txParams.T_t_tx=rfParams.T_t*(1+txParams.T_t_err/100);
%% 
% determine suitable oversampling/downsampling factor

[M,N]=rat(rfParams.f_s*txParams.T_t_tx);

disp('oversampling/downsampling factor: ');
disp(['M = ', num2str(M),', N = ', num2str(N)]);
%% 
% upsample the message signals

ms4up=upsample(ms4,M);
mc4up=upsample(mc4,M);
%% 
% generate pulse-shaped signals

psf=srrc(rfParams.SRRCLength,rfParams.SRRCrolloff,M,0);
xs=conv(psf',ms4up)';
xc=conv(psf',mc4up)';

fig=1; % Plot pulse-shaping filter response
if plotsONOFF(fig)==1
    figure(fig)
    freqz(psf,512,0:(pi/8)/512:pi/8), zoom xon
    subplot(2,1,1),title('Pulse-Shaping Filter Frequency Response')
end

%% 
% create time vectors

% create time vector, t
Ts=txParams.T_t_tx/M;
t=0:Ts:Ts*(length(ms4up)-1);

% t range for plot
% tRange=1:length(t); % all time
t1=floor((length(frameParams.preamble)+floor(frameParams.userDataLength/4)+10)*4*M);
tRange=t1:t1+25*M;  % t range for plot

%%
fig=fig+1; % Plot 1/2 of Message as 2-PAM symbols
if plotsONOFF(fig)==1
    figure(fig)
    square_s=repelem(ms4,M);                % oversampled ms4 as square pulses
    plotspec4(square_s,Ts,t,tRange), zoom xon
    axis([-1/(Ts*M) 1/(Ts*M),ylim])
    subplot(2,1,1), axis([xlim,1.1*ylim])
    subplot(2,1,1), title('1/2 of Message as 2-PAM Symbols')
end    
    
fig=fig+1; % Plot 1/2 of Message as 2-PAM symbols
if plotsONOFF(fig)==1
    figure(fig)
    square_c=repelem(mc4,M);                % oversampled mc4 as square pulses
    plotspec4(square_c,Ts,t,tRange)
    axis([-1/(Ts*M) 1/(Ts*M),ylim])
    subplot(2,1,1), axis([xlim,1.1*ylim]), hold off
    subplot(2,1,1), title('Other 1/2 of Message as 2-PAM symbols')
end

t=0:Ts:Ts*(length(xs)-1);
tRange=tRange+floor(0.5*rfParams.SRRCLength*M);
fig=fig+1; % Plot Pulse-Shaped 1/2 Message Signal
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(xs,Ts,t,tRange), zoom xon, axis([-1/(Ts*M) 1/(Ts*M),ylim])
    subplot(2,1,1), title('Pulse-Shaped 1/2 Message Signal')
end
%% 
% CREATE CARRIERS

fc=txParams.f_if_tx;

% phase offset
if phaseErrorOn == 1
    th=-1;
else    
    th=0;
end

% phase noise 
if phaseNoiseON == 1
    p_noise=cumsum(randn(size(xs))*sqrt(txParams.phaseNoiseVariance/N));
end

% ADD PHASE NOISE AND/OR PHASE OFFSET TO CARRIERS
if phaseNoiseON == 1
    % generate phase noise process
    p_noise=cumsum(randn(size(xs))*sqrt(txParams.phaseNoiseVariance/N));
    % create carrier with phase noise
    cc=cos(2*pi*fc*t+th+p_noise');
    cs=sin(2*pi*fc*t+th+p_noise');

    fig=fig+1;  % plot phase noise
    if plotsONOFF(fig)==1
        figure(fig)
        plotspec4(p_noise,Ts,t,tRange), zoom xon
        subplot(2,1,1), title('Phase Noise Signal Added to Carriers')
    end
else % no phase noise
    cc=cos(2*pi*fc*t+th);
    cs=sin(2*pi*fc*t+th);
    fig=fig+1;  % skip phase noise plot
end

fig=fig+1; % Plot Transmitter Cosine Carrier Signal
tRange=floor(length(t)/2):floor(length(t)/2)+M;
fLim=min([2*rfParams.f_if,0.5*M/txParams.T_t_tx]);
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(cc,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1), title('Transmitter Cosine Carrier Signal')
end
%% 
% mix signal to RF (well, technically IF)

x_rf=xc.*cc'-xs.*cs'; % Quadrature Amplitude Modulation (QAM)

tRange=floor(length(t)/2):floor(length(t)/2)+25*M;
fig=fig+1; % Plot Transmitter Cosine Carrier Signal
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(x_rf,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1), title('Transmitter QAM Signal'),saveylim=ylim;
end
%% 
%CHANNEL IMPAIRMENTS
% add multipath

mpf=upsample(chanParams.multipath,M);    % multipath impulse response
r=filter(mpf,1,x_rf);    % output of channel

fig=fig+1; % Plot multipath impulse response
if plotsONOFF(fig)==1
    figure(fig)
    freqz(mpf), zoom xon
    subplot(2,1,1),title('Multipath Impulse Response')
    subplot(2,1,2),plot(mpf,'-o')
end

fig=fig+1; % Plot Signal with Multipath
if plotsONOFF(fig)==1
    figure(fig)
    %tRange=tRange+floor(0.5*numel(chanParams.multipath)*M);
    plotspec4(r,Ts,t,tRange), zoom xon, axis([-fLim,fLim,ylim])
    subplot(2,1,1), title('Signal with Multipath'), axis([xlim,saveylim])
end
%% 
% add noise

sp=pow(r);
np=sp/10^(chanParams.SNR/10);
n=sqrt(np)*randn(1,length(r));
r=r'+n;

fig=fig+1; % Plot Signal Plus Noise
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(r,Ts,t,tRange), zoom xon
    axis([rfParams.f_if-2/rfParams.T_t rfParams.f_if+2/rfParams.T_t,ylim])
    subplot(2,1,1), title('Signal Plus Noise'), axis([xlim,saveylim])
end

%% 
% perform downsampling by N, normalization

r=r(1:N:end);
r=r(:)*sqrt(length(r))/norm(r);

Ts=N*txParams.T_t_tx/M;
t=0:Ts:Ts*(length(r)-1);
tRange=floor(length(t)/2):floor(length(t)/2)+25*floor(M/N); % middle
% tRange=1:25*floor(M/N);                                   % beginning

fig=fig+1; % Plot Transmitted QAM Signal
if plotsONOFF(fig)==1
    figure(fig)
    plotspec4(r,Ts,t,tRange), zoom xon%, axis([-fLim,fLim,ylim])
    axis([rfParams.f_if-2/rfParams.T_t rfParams.f_if+2/rfParams.T_t,ylim])
    subplot(2,1,1), title('Received and Downsampled QAM Signal')
end

clearvars -except rfParams frameParams r mc4 ms4