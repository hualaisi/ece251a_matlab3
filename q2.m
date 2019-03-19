close all;clc;clear

% filter parameter
num = [1, -0.9, 0.81]; % numerator 
den = [1, -2.76, 3.809, -2.654, 0.924]; % denominator 
variance = 1; % variance
samples = [128, 256];
poverlap = [0.5, 0.75]; % overlap percentage
nsegments = [32, 64, 128]; % # segments 

% Gaussian white noise with variance 1
%% hamming window
for i = 1:length(samples)
    psd = psdTrue(samples(i), num, den);
    figure
    count = 1;
    for j = 1:length(nsegments)
        for k = 1:length(poverlap)
            temp = 0;
            iteration = 100;
            for t = 1:iteration
                [pxx,w] = pwelch(gwn(samples(i), num, den), nsegments(j), nsegments(j)*poverlap(k));
                temp = temp + pxx;
            end
            temp = temp/iteration;
            
            subplot(3,2,count);
            plot(w, 10*log10(abs(temp)))
            xlim([0,pi])
            hold on
            
            freq = 0:samples(i)/2-1;
            plot(2*freq/samples(i)*pi,10*log10(psd), '--')
            xlim([0,pi])
            hold off
            legend('Welch Method PSD','True PSD')
            grid on
            title(['Hamming; Sample: ',num2str(samples(i)),', Segments: ',num2str(nsegments(j)),', Overlap: ',num2str(poverlap(k))])
            count = count + 1;
        end
    end
end

%% rectangular window
for i = 1:length(samples)
    psd = psdTrue(samples(i), num, den);
    figure
    count = 1;
    for j = 1:length(nsegments)
        for k = 1:length(poverlap)
            temp = 0;
            for t = 1:100
                HH = rectwin(nsegments(j));
                [pxx,w] = pwelch(gwn(samples(i), num, den), HH, nsegments(j)*poverlap(k));
                temp = temp + pxx;
            end
            temp = temp/1000;
            
            subplot(3,2,count);
            plot(w, 20*log10(abs(temp)))
            xlim([0,pi])
            hold on
            
            freq = 0:samples(i)/2-1;
            plot(2*freq/samples(i)*pi,10*log10(psd), '--')
            xlim([0,pi])
            hold off
            legend('Welch Method PSD','True PSD')
            grid on
            title(['Rectangular: Sample: ',num2str(samples(i)),', Segments: ',num2str(nsegments(j)),', Overlap: ',num2str(poverlap(k))])
            count = count + 1;
        end
    end
end

%% functions
% fucntion to calculate psd
function psd = psdTrue(sample, num, den) 
    % get true autpcorrelation
    xn_impulse = [1, zeros(1, sample-1)];
    yn_impulse = filter(num, den, xn_impulse);
    
    % get true psd
    psd = fft(yn_impulse, sample);
    psd = abs(psd).*abs(psd);
    psd = abs(psd(1:sample/2));
end

% fucntion to get autocorrelation of Gaussian white noise with variance 1
function yn_est = gwn(sample, num, den)
    variance = 1;
    xn_est = sqrt(variance)*wgn(1, sample, 0); % xn estimate
    yn_est = filter(num, den, xn_est); % yn estimate
end