function [averagedestimate] = avgest(signal,M)
%{AVGEST Averages the PSD for a less noisy estimate by splitting up the signal into M segments and then calculates the ACF for each part and averages the PSD.
%
%   output: averagedestimate    the resulting averagedestimate
%
%   input:  signal              The signal which PSD will be averaged
%           M                   Determines how many parts the periodogram
%                               will split into
%}

N = max(size(signal));


if (M > sqrt(N))
    disp('Warning: Distortions may occur due to high M')
end

n = N/M;

avgacf = zeros(n,n);

for k = 1:M
    avgacf(1:n,k) = bartlett(signal((k-1)*n+1:k*n));
end

avgacf = avgacf';

%quickfix for scaling issues that only seem to appear when M<sqrt(N)
if (M < sqrt(N))
avgacf = avgacf*(N/M^2); 
end

averagedestimate = abs(fft(mean(avgacf)));
end

