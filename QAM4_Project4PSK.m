%% INFO and SETUP
% r = input signal
% Ts = sampling rate
% M = Oversampling factor
% N = Downsampling factor
clear all; 
close all;

SDR_Proj_4QAM_Tx_R2;
M = rfParams.f_s*rfParams.T_t;
fs = rfParams.f_s;
Ts = 1/fs;
t = Ts:Ts:length(r)*Ts;
fc = rfParams.f_if;
tRange = floor(length(t)/2):floor(length(t)/2)+25*floor(M); 
 
fig = 0;
fig = fig+1;
figure(fig)
plotspec4(r,Ts,t,tRange), zoom xon, axis([rfParams.f_if-2/rfParams.T_t rfParams.f_if+2/rfParams.T_t,ylim]);
%plotspec(r,Ts);
title('Received Signal');
%% Costas Loop
OverSamp = M;
fl=200;
fbe=2*[0 fc*Ts-1/(16*OverSamp) fc*Ts+1/(16*OverSamp) .5];
amps=[1 1 0 0];
b=firls(fl,fbe,amps);

recSig = r';

fig = fig+1;
figure(fig)
plotspec(b,Ts);
cmu1 = .001;
cmu2 = .0004;

theta=zeros(1,length(recSig));
theta2=zeros(1,length(recSig));

v1=zeros(1,length(recSig));
v2=zeros(1,length(recSig));
v3=zeros(1,length(recSig));
v4=zeros(1,length(recSig));

v1=filter(b,1,recSig.*sin(2*pi*fc*t));
v2=filter(b,1,recSig.*sin(2*pi*fc*t+pi/4));
v3=filter(b,1,recSig.*cos(2*pi*fc*t));
v4=filter(b,1,recSig.*cos(2*pi*fc*t+pi/4));

Bdelay=(length(b)-1)/2;
z1=zeros(1);
z2=zeros(1);
z3=zeros(1);
z4=zeros(1);
for k=1:length(recSig)-1-length(b)/2
  cosThetaK=cos(theta(k));
  sinThetaK=sin(theta(k));
    z1=v1(k+Bdelay)*cosThetaK+v3(k+Bdelay)*sinThetaK;
    z2=v2(k+Bdelay)*cosThetaK+v4(k+Bdelay)*sinThetaK;
    z3=v3(k+Bdelay)*cosThetaK-v1(k+Bdelay)*sinThetaK;
    z4=v4(k+Bdelay)*cosThetaK-v2(k+Bdelay)*sinThetaK;
    theta(k+1)=theta(k)+cmu1*(z2*z4*z3*z1);

    cosTheta12=cos(theta(k)+theta2(k));
    sinTheta12=sin(theta(k)+theta2(k));
    z1=v1(k+Bdelay)*cosTheta12+v3(k+Bdelay)*sinTheta12;
    z2=v2(k+Bdelay)*cosTheta12+v4(k+Bdelay)*sinTheta12;
    z3=v3(k+Bdelay)*cosTheta12-v1(k+Bdelay)*sinTheta12;
    z4=v4(k+Bdelay)*cosTheta12-v2(k+Bdelay)*sinTheta12;
    theta2(k+1)=theta2(k)+cmu2*(z4*z2*z3*z1);
end

fig = fig+1;
figure(fig)
subplot(2,1,1),plot(theta);
subplot(2,1,2),plot(theta2);

downCSig=zeros(1,length(recSig));
downCSig=2*recSig.*exp(-j*(2*pi*fc*t+theta+theta2));
downSig = filter(b,1,downCSig);

fig = fig+1;
figure(fig)
plotspec4(downSig,Ts,t,tRange), zoom xon, axis([-2/rfParams.T_t 2/rfParams.T_t,ylim]);
title('LPF Costas Output');

%% Timing Recovery
ApproxPeriods = rfParams.SRRCLength;
beta = rfParams.SRRCrolloff;
Time_mu1 = 0.0010;
Time_mu2 = 0.000025;

stuff=(ApproxPeriods*OverSamp+1):OverSamp:(length(downSig)-ApproxPeriods*OverSamp);

tau=zeros(1,length(stuff));
tau2=zeros(1,length(stuff));
timRecSig=zeros(1,length(stuff));
clear stuff;

i=1;
delta=0.00001;

timRecPDelta=2*ApproxPeriods*OverSamp+1;
timRecMDelta=2*ApproxPeriods*OverSamp+1;
timRecPDelta2=2*ApproxPeriods*OverSamp+1;
timRecMDelta2=2*ApproxPeriods*OverSamp+1;

for k=(ApproxPeriods*OverSamp+1):OverSamp:(length(downSig)-ApproxPeriods*OverSamp)
  data=downSig((k-ApproxPeriods*OverSamp):(k+ApproxPeriods*OverSamp));
  
  timRecSig(i)=srrc(ApproxPeriods,beta,OverSamp,tau(i)+tau2(i))*data.';
  
  tempSig=srrc(ApproxPeriods,beta,OverSamp,tau(i))*data.';
  
  timRecPDelta=srrc(ApproxPeriods,beta,OverSamp,tau(i)+delta)*data.';
  
  timRecMDelta=srrc(ApproxPeriods,beta,OverSamp,tau(i)-delta)*data.';
  
  timRecPDelta2=srrc(ApproxPeriods,beta,OverSamp,tau(i)+tau2(i)+delta)*data.';
  
  timRecMDelta2=srrc(ApproxPeriods,beta,OverSamp,tau(i)+tau2(i)-delta)*data.';

  dXdTau2=(abs(timRecPDelta2)-abs(timRecMDelta2))/delta;
  tau2(i+1)=tau2(i)-Time_mu2*abs(timRecSig(i))^3*dXdTau2;
  
  dXdTau=(abs(timRecPDelta)-abs(timRecMDelta))/delta;
  tau(i+1)=tau(i)-Time_mu1*abs(tempSig)^3*dXdTau;  
  
  i=i+1;
end

fig = fig+1;
figure(fig);
subplot(2,1,1),plot(timRecSig,'.b');
subplot(2,1,2),plot(tau+tau2);
title('Timing Recovery');
%% LMS Equalizer
%% Quantize
decOut_re=sign(real(timRecSig));
decOut_im=-sign(imag(timRecSig));
%decOut=decOut(length(training)+1+delta:length(training)+length(msg1sc)+delta);
%% Frame Sync
preamble = frameParams.preamble;
datalen = frameParams.userDataLength*8;

header = zeros(1,8*length(preamble));
for k=0:length(preamble)-1                  % code characters in 8-bit ascii
  bc=dec2base(double(preamble(k+1)),2,8);
  for i=0:7
    header(8*k+1+i)=2*bin2dec(bc(i+1))-1;
  end  
end

y1=xcorr(header , decOut_re);
[mymax1, ind1 ]=max(abs(y1));
datastart1=mod(length(decOut_re)-ind1 , length(decOut_re) )+1;

y2=xcorr(header , decOut_im);
[mymax2, ind2 ]=max(abs(y2));
datastart2=mod(length(decOut_im)-ind2 , length(decOut_im) )+1;

fig=fig+1;
figure(fig)
subplot(2,1,1),plot(y1);
subplot(2,1,2),plot(y2);
title('Correlations');

%% Decode
offset = 1027+length(header);
m_out1 = zeros(1,length(decOut_im-offset));
m_out2 = zeros(1,length(decOut_re-offset));
k = 1;

for i = offset:4:length(decOut_im)-4
    b0 = decOut_im(i);
    b1 = decOut_im(i+1);
    b2 = decOut_im(i+2);
    b3 = decOut_im(i+3);
    m_out1(k) = qpsk_decode(b0,b1);
    m_out1(k+1)= qpsk_decode(b2,b3);
    k = k+2;
end
k=1;
for i = offset:4:length(decOut_re)-4
    b0 = decOut_re(i);
    b1 = decOut_re(i+1);
    b2 = decOut_re(i+2);
    b3 = decOut_re(i+3);
    m_out2(k) = qpsk_decode(b0,b1);
    m_out2(k+1)= qpsk_decode(b2,b3);
    k = k+2;
end

message_out1 = pam2letters(m_out1);
message_out2 = pam2letters(m_out2);

m_outf = char(zeros(1,2*length(message_out2)));
k=1;
for i = 1:length(message_out2)-1
    m_outf(k)=message_out2(i);
    m_outf(k+1)=message_out1(i);
    k=k+2;
end

message=m_outf


