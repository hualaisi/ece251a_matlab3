close all;clc;clear

% filter parameter
num = [1, -0.9, 0.81]; % numerator 
den = [1, -2.76, 3.809, -2.654, 0.924]; % denominator 
variance = 1; % variance
noiseLength = 1024; % length gaussian white noise 
samples = [64, 128, 256, 512, 1024];


%% part (a)
% get true autpcorrelation
xn_impulse = [1, zeros(1, 1023)];
yn_impulse = filter(num, den, xn_impulse);

figure
plot(autocorr(yn_impulse, noiseLength-1))
grid on
xlim([0, noiseLength-1])
xlabel('Index m');
title('True Autocorrelation Function')


%% part (b)
xn_est = sqrt(variance)*wgn(1, noiseLength, 0); % xn estimate
yn_est = filter(num, den, xn_est); % yn estimate

% get autocorrelation 

for n = 1:length(samples)
    figure
    subplot(2,1,1)
    plot(autocorr(yn_est, samples(n)-1))
    grid on
    xlim([0, noiseLength-1])
    xlabel('Index m')
    title(['Estimated Autocorrelation with ',num2str(samples(n)),' Samples'])
    
    subplot(2,1,2)
    plot(autocorr(yn_impulse, noiseLength-1))
    grid on
    xlim([0, noiseLength-1])
    xlabel('Index m');
    title('True Autocorrelation Function')
end


%% part (c)
Fs = 1024;
N = length(yn_impulse);
xdft = fft(yn_impulse);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(yn_impulse):Fs/2;

figure
plot(freq/Fs*2*pi,psdx)
grid on
xlabel('Frequency (Rad/s)')
ylabel('Power Spectral Density (W/Hz)')
title('True Power Spectrum')


%% part (d)
xn_est_d = sqrt(variance)*wgn(1, noiseLength, 0); % xn estimate
yn_est_d = filter(num, den, xn_est_d); % yn estimate

% get autocorrelation 
for n = 1:length(samples)
    Fs = 1024;
    N = samples(n);
    xdft_d = fft(yn_est_d(1:samples(n)));
    xdft_d = xdft_d(1:N/2+1);
    psdx_d = (1/(Fs*N)) * abs(xdft_d).^2;
    psdx_d(2:end-1) = 2*psdx_d(2:end-1);
    freq_d = 0:Fs/length(yn_est_d(1:samples(n))):Fs/2;

    figure
    subplot(2,1,1)
    plot(freq_d/Fs*2*pi,psdx_d/Fs)
    grid on
    xlabel('Frequency (Rad/s)')
    ylabel('Power Spectral Density (W/Hz)')
    title(['Periodogram of Estimated Autocorrelation with ',num2str(samples(n)),' Samples'])
    
    subplot(2,1,2)
    plot(freq/Fs*2*pi,psdx)
    grid on
    xlabel('Frequency (Rad/s)')
    ylabel('Power Spectral Density (W/Hz)')
    title('Periodogram of True Autocorrelation Function')
end


%% part (e)
iteration = 1000; % numbers of iterations

for n = 1:length(samples)
    sample_array = []; % create an empty array for each samples
    
    for i = 1:iteration
        xn_est_e = sqrt(variance)*wgn(1, noiseLength, 0); % xn estimate
        yn_est_e = filter(num, den, xn_est_e); % yn estimate

        Fs = 1024;
        N = samples(n);
        xdft_e = fft(yn_est_e(1:samples(n)));
        xdft_e = xdft_e(1:N/2+1);
        psdx_e = (1/(Fs*N)) * abs(xdft_e).^2;
        psdx_e(2:end-1) = 2*psdx_e(2:end-1);
        sample_array = [sample_array; psdx_e];
    end
    
    freq_e = 0:Fs/length(yn_est_e(1:samples(n))):Fs/2;
        
    % plot of mean
    figure
    plot(freq_e/Fs*2*pi,mean(sample_array)/Fs)
    hold on
    plot(freq/Fs*2*pi,psdx)
    hold off
    legend('Estimated', 'True Autocorrelation');
    
    grid on
    xlabel('Frequency (Rad/s)')
    ylabel('Power Spectral Density (W/Hz)')
    title(['Mean of Periodogram of Estimated Autocorrelation with ',num2str(samples(n)),' Samples'])
    
    % plot of variance
    figure
    plot(freq_e/Fs*2*pi,var(sample_array))
    grid on
    xlabel('Frequency (Rad/s)')
    title(['Variance of Periodogram of Estimated Autocorrelation with ',num2str(samples(n)),' Samples'])
    
end
