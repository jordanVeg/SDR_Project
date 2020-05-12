%%INFO and SETUP
% r = input signal
% Ts = sampling rate
% M = Oversampling factor
% N = Downsampling factor
clear all; 
close all;

SDR_Proj_4PAM_Tx;
 
fig = 0;
fig = fig+1;
figure(fig)
plotspec4(r,Ts,t,tRange), zoom xon, axis([1800,2200,ylim]);
title('Received Signal');
M = rfParams.f_s*rfParams.T_t;
% fs = rfParams.f_s*M;
% Ts = 1/fs;
%% Downconvert and Carrier Recovery
fc = 2000;
offset=0.0012;

dmod=recover_sig(r,t,rfParams.f_s,Ts,offset);
% fig = fig+1;
% figure(fig)
% plotspec(dmod,Ts)
% title('Recovered Carrier Signal');

carrier = cos(2*pi*fc*t);
dmod = carrier;
% fig = fig+1;
% figure(fig)
% plot(dmod - carrier);
% title('Difference of Carrier Signals');

rec = r.*dmod;

fig = fig+1;
figure(fig)
%plotspec(rec,Ts);
plotspec4(rec,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
title('Downconverted Signal');

% fig = fig+1;
% figure(fig)
% subplot(2,1,1), plot(th1);
% subplot(2,1,2), plot(th2);
% title('Freq-Phase Tracking');
%% Match Filter
% LPF
bandstart = 0.5;
bandend = bandstart + bandstart*0.12;
freqs=[0 0.4 0.45 1];
amps=[1 1 0 0];
f0 = 100;
b=firpm(f0,freqs,amps);      
lpf_out = filter(b,1,rec);  
fig = fig+1;
figure(fig);
plotspec4(lpf_out,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
title('Lowpass Filter Signal');
srrc_length = rfParams.SRRCLength;

src = srrc(srrc_length,rfParams.SRRCrolloff,M,0);


matched = conv(rec,src);
fig = fig+1;
figure(fig);
plotspec4(matched,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
title('Matched Filter Signal');


%% Interpolater Downsampler and Timing Recovery
 M_r=round(rfParams.f_s*rfParams.T_t);
offset = 1 + 9*M_r;
[rdown,t_error] = recover_timing(matched(offset:end),srrc_length,M,Ts);
% fig = fig+1;
% figure(fig)
% %plot([1:length(tr)],tr,'.');
% plotspec4(tr,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
% title('Timing Recovered Signal')


% % approximates the upsample rate to nearest integer
% M1=round(rfParams.f_s*rfParams.T_t);
% offset=9*M1;  % offset to get to the first symbol of the message (9 is for the 4-PAM Tx signal)
% sd=zeros(1,length(tr)); % z is the baseband signal with sampling rate of f_s; sd is vector to save the samples T_t apart
% indexes=1+offset:rfParams.f_s*rfParams.T_t:length(tr);  % exact evenly spaced fractional indices T_t spaced 
% indexes=round(indexes); % indexes rounded to nearest integer
% sd(indexes)=tr(indexes);  % save the values of z for the almost evenly spaced time samples T_t apart 
% rdown=tr(indexes); % remove the zeros in between the T_t spaced samples; z is downsampled by the unrounded M value and placed in z2

fig = fig+1;
figure(fig);
subplot(2,1,1), plot([1:length(rdown)],rdown,'.');
subplot(2,1,2), plot([1:length(t_error)],t_error);
title('Downsampled Signal');


%% Equalizer
% len = length(rdown);
% n=6;                        % equalizer filter length
% f=zeros(1,n); f(floor(n/2))=1; f=f';  % center spike initialization
% mu=0.0005;%0.003 - 0.012      % algorithm stepsize
% em=zeros(1,len);               % error matrix
% fm=zeros(n,len);               % filter matrix
% for i=n+1:len                  % iterate
%   rr=r(i:-1:i-n+1)';         % vector of received signal
%   e=quantalph(f'*rr, [-3 -1 1 3])-f'*rr;       % calculate error
%   em(i)=e;                   % save error value
%   f=f+mu*e*rr;               % update equalizer coefficients
%   fm(:,i)=f;                 % save filter coefficients
% end
% rdown = filter(f,1,rdown);
% fig = fig+1;
% figure(fig);
% plot([1:length(rdown)],rdown,'.');
% title('Equalized Signal');
%% AGC
n = length(rdown);
ds=sqrt(20);                       % desired power of output
mu=0.00018;                          % algorithm stepsize
lenavg=100;                         % length over which to average
a=zeros(n,1); a(1)=2;             % initialize AGC parameter
sp=zeros(n,1);                      % initialize outputs
avec=zeros(1,lenavg);              % vector to store terms for averaging
for k=1:n-1
  sp(k)=a(k)*rdown(k);                  % normalize by a(k)
  avec=[sign(a(k))*(sp(k)^2-ds),avec(1:lenavg-1)];  % incorporate new update into avec
  a(k+1)=a(k)-mu*mean(avec);       % average adaptive update of a(k)
end

fig = fig+1;
figure(fig)
subplot(2,1,1),plot([1:length(sp)],sp,'.');
subplot(2,1,2),plot(a);
title('Downsampled with Gain');
%% Decision Device
alphabet = [-3 -1 1 3];
dec = quantalph(sp,alphabet)';

fig = fig+1;
figure(fig)
plot([1:length(dec)],dec,'.');
title('Quantized Signal');
%% Decoder
%temp = circshift(dec,[-1,0]);
temp = zeros(size(dec));
for i = 1:length(dec)-1
    temp(i+1)= dec(i);
end
temp(1)= dec(end);
message_out = pam2letters(temp);