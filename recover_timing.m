function [recovered_sig,error]=recover_timing(x,l,m,Ts)

n = length(x(1:round(m):end));
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.1;                            % algorithm stepsize
delta=0.1;                    % time for derivative
while tnow<length(x)-2*l*m            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  tau=tau+mu*dx*xs(i);              % alg update (energy)
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end
recovered_sig = xs(1:i-2);
error = tausave(1:i-2);
end
