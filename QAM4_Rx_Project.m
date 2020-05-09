%%INFO and SETUP
% r = input signal
% Ts = sampling rate
% M = Oversampling factor
% N = Downsampling factor
fig = 0;
fig = fig+1;
figure(fig)
%plotspec(r,Ts);
plotspec4(r,Ts,t,tRange), zoom xon, axis([1800,2200,ylim]);
title('Received Signal');
%% Downconvert
fc = 2000;
dmod = cos(2*pi*fc*t);
rec = r.*dmod;

fig = fig+1;
figure(fig)
%plotspec(rec,Ts);
title('Downconverted Signal');
plotspec4(rec,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
%% Match Filter
% LPF
bandstart = 0.5;
bandend = bandstart + bandstart*0.12;
freqs=[0 bandstart bandend 1];
amps=[1 1 0 0];
f0 = 100;
b=firpm(f0,freqs,amps);      
%lpf_out = filter(b,1,rec);   
srrc_length = rfParams.SRRCLength;

src = srrc(srrc_length,rfParams.SRRCrolloff,rfParams.f_s*rfParams.T_t,0);

%normalizing filter
%nf = filter(fliplr(src)/pow(src)*(N/M),1,lpf_out);
%nf = lpf_out;

nf = conv(rec,src);
fig = fig+1;
figure(fig);
plotspec4(nf,Ts,t,tRange), zoom xon, axis([-4000,4000,ylim]);
title('LPF Signal');
% fig = fig+1;
% figure(fig)
% plotspec(b,Ts)
%% Interpolater Downsampler

  % approximates the upsample rate to nearest integer
  M1=round(rfParams.f_s*rfParams.T_t);
offset=9*M1;  % offset to get to the first symbol of the message (9 is for the 4-PAM Tx signal)
sd=zeros(1,length(nf)); % z is the baseband signal with sampling rate of f_s; sd is vector to save the samples T_t apart
indexes=1+offset:rfParams.f_s*rfParams.T_t:length(nf);  % exact evenly spaced fractional indices T_t spaced 
indexes=round(indexes); % indexes rounded to nearest integer
sd(indexes)=nf(indexes);  % save the values of z for the almost evenly spaced time samples T_t apart 
rdown=nf(indexes); % remove the zeros in between the T_t spaced samples; z is downsampled by the unrounded M value and placed in z2

fig = fig+1;
figure(fig);
plot([1:length(rdown)],rdown,'.');
title('Downsampled Signal');

% fig = fig+1;
% figure(fig);
% plotspec(rdown,1/Ts)
% title('Downsampled Pulse Shape');
% x = rdown;
% n = length(x);
% l = (srrc_length);
% tnow=l*M+1; tau=0; xs=zeros(1,n);   % initialize variables
% tausave=zeros(1,n); tausave(1)=tau; i=0;
% mu=0.005;                            % algorithm stepsize
% delta=0.1;                          % time for derivative
% while tnow<length(x)-l*M            % run iteration
%   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
%   dx=x_deltap-x_deltam;             % numerical derivative
%   tau=tau+mu*dx*xs(i);              % alg update (energy)
%   tnow=tnow+m; tausave(i)=tau;      % save for plotting
% end
% rdown = xs;
% 
% fig = fig+1;
% figure(fig);
% subplot(2,1,1)
% plot(rdown(1:i-2),'.');
% subplot(2,1,2), plot(tausave(1:i-2))   
% title('DDClockRec Signal');

%AGC
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
plot([1:length(sp)],sp,'.');
title('Downsampled with Gain');

%% Equalizer
% Not needed for perfect channel
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