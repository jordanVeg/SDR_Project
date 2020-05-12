% plotspec4(x,Ts,t,tRange) plots the spectrum of the signal x
% Ts = time (in seconds) between adjacent samples in x
% t = vector of sample time values
% tRange = Range of indexes to zoom into for the time domain plot of x
function plotspec4(x,Ts,t,tRange)
N=length(x);                               % length of the signal x
%t=Ts*(1:N);                                % define a time vector
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);       % frequency vector
fx=fft(x(1:N));                            % do DFT/FFT
fxs=fftshift(fx);                          % shift it for plotting
subplot(2,1,1), plot(t(tRange),x(tRange))  % plot the waveform
xlabel('seconds'); ylabel('amplitude')     % label the axes
subplot(2,1,2), plot(ssf,abs(fxs))         % plot magnitude spectrum
xlabel('frequency'); ylabel('magnitude')   % label the axes