function carrier = recover_sig(r,t,fs,Ts,offset)
% r: input signal to be recovered
% t: input time vector
%% PreProcess PLL
q=r.^2;                           % square nonlinearity
fl=500; ff=[0 .899 .90 .99 .999 1]; % BPF center frequency at .4
%ff = [0 850 950 1020 1150 fs/2]*2/fs;
fa=[0 0 1 1 0 0];                 % which is twice f_0
h=firpm(fl,ff,fa);                % BPF design via firpm
rp=filter(h,1,q);      % filter gives preprocessed r

%% PLL
mu= offset;                              % algorithm stepsize
a=[-1 1]; lena=length(a)-1;              % autoregressive coefficients
b=20*[-2 2-mu]; lenb=length(b);             % moving average coefficients
xvec=zeros(lena,1); evec=zeros(lenb,1);  % initial states in filter
f_Tx=2000.0;                               % assumed freq. at receiver
theta=zeros(1,length(t)); theta(1)=0;    % initialize vector for estimates
for k=1:length(t)-1                      % e contains past fl+1 inputs
  e=rp(k)*sin(4*pi*f_Tx*t(k)+2*theta(k)); % input to filter
  evec=[e;evec(1:lenb-1)];               % past values of inputs
  x=-a(2:lena+1)*xvec+b*evec;            % output of filter
  xvec=[x;xvec(1:lena-1,1)];             % past values of outputs
  theta(k+1)=theta(k)+mu*x;              % algorithm update
end
carrier = cos(2*pi*f_Tx*t+theta(end));
figure(10)
plotspec(carrier,Ts)
end
