%{
Study 1 for the course TSDT14 at Linkoping University
Estimates the ACF for a low grade LP filter as well as an ideal LP filter.
Then plots them with the theoretical values to give a comparison

Author: Hoang Tran, hoatr725 (2020)
%}

clear;
clc;
close all;

%Noise construction and misc.
N = 2^16;
x = randn(N,1); %noise
ACFplotSpace = -N/2+1:1:N/2;

%noise plot
%figure
%plot(0:1:N, x)

%% low grade filter 3 elements averaging filter
B = ones(3,1)/3; %filter kernel
lowGradeOutput = filter(B, 1, x); %filtered with B
lgACF = bartlett(lowGradeOutput);
lgPSD = fft(lgACF);
PSDplotSpace = linspace(0,1,length(lgPSD));
tlgPSD = ((1/3 * (1+2*cos(2*pi*PSDplotSpace))).^2)'; %theoretical PSD form from modulo of frequency response squared
tlgACF = circshift(fftshift(ifft(tlgPSD)),-1); %theoretical ACF shifted


figure  %figure for low degree filter
subplot(2,2,1)
plot (ACFplotSpace,lgACF)
title('Estimated ACF of low grade LP filter');
xlabel('k');
ylabel('ry[k]');
subplot(2,2,2)
plot (PSDplotSpace, abs(lgPSD))
title('Estimated PSD of low grade LP filter');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot (2,2,4)
plot(PSDplotSpace, abs(tlgPSD))
title('Theoretical PSD of low grade LP filter');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot (2,2,3)
plot(ACFplotSpace, tlgACF)
title('Theoretical ACF of low grade LP filter');
xlabel('k');
ylabel('ry[k]');

%stem plot
figure
subplot (2,1,1)
stem(ACFplotSpace, lgACF); xlim([-20,20]);
title('Stem plot of estimated ACF of low grade LP filter');
xlabel('k');
ylabel('ry[k]');

subplot (2,1,2)
stem(ACFplotSpace, tlgACF); xlim([-20,20]);
title('Stem plot of theoretical ACF of low grade LP filter');
xlabel('k');
ylabel('ry[k]');

%% Ideal filter with cutoff frequency
cutOffFrequency = 0.1;
[IFB,IFA] = butter(30, cutOffFrequency*2);
idealOutput = filter(IFB, IFA, x);
idealACF = bartlett (idealOutput);
idealPSD = fft(idealACF);

%theoretical ideal filter PSD
tidealPSD = zeros(N,1);
tidealPSD(1:floor(cutOffFrequency*N)) = 1;
tidealPSD(ceil((1-cutOffFrequency)*N+1):N) = 1;
tidealACF = circshift(fftshift(ifft(tidealPSD)),-1); %theoretical ACF shifted

figure
subplot(2,2,1)
plot (ACFplotSpace,idealACF)
title('Estimated ACF of high grade LP filter');
xlabel('k');
ylabel('ry[k]');
subplot(2,2,2)
plot (PSDplotSpace, abs(idealPSD))
title('Estimated PSD of high grade LP filter');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot (2,2,4)
plot(PSDplotSpace, abs(tidealPSD))
title('Theoretical PSD of high grade LP filter');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot (2,2,3)
plot(ACFplotSpace, tidealACF)
title('Theoretical ACF of high grade LP filter');
xlabel('k');
ylabel('ry[k]');

%stem plot
figure
subplot (2,1,1)
stem(ACFplotSpace, idealACF); xlim([-20,20]);
title('Stem plot of estimated ACF of high grade LP filter');
xlabel('k');
ylabel('ry[k]');
subplot (2,1,2)
stem(ACFplotSpace, tidealACF); xlim([-20,20]);
title('Stem plot of theoretical ACF of high grade LP filter');
xlabel('k');
ylabel('ry[k]');
