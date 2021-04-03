%{
Study 3 for the course TSDT14 at Linkoping University
Study of Non-LTI-systems
Author: Hoang Tran, hoatr725 (2020)
%}

clear;
clc;
close all;
%%

%Noise construction and misc.
N = 2^16;
x = randn(N,1); %noise
PSDplotSpace = linspace(0,1,N);

%filtering
cutOffFrequency = 0.1;
[IFB,IFA] = butter(30, cutOffFrequency*2);
FilterOutput = filter(IFB, IFA, x);

%ideal filter
tidealPSD = zeros(N,1);
tidealPSD(1:floor(cutOffFrequency*N)) = 1;
tidealPSD(ceil((1-cutOffFrequency)*N+1):N) = 1;

figure
subplot(2,1,1)
plot(PSDplotSpace, tidealPSD)
title('Theoretical PSD of output from high order LP filter');
xlabel('\theta');
ylabel('Rx(\theta)');
subplot(2,1,2)
histogram(FilterOutput)
title('Histogram of the filter output');
xlabel('Amplitude');
ylabel('Occurrences');
%% theoretical and estimated PSD of squarer output
%theoretical output
tPSDsq =  0.4*(tripuls(PSDplotSpace/(4*cutOffFrequency))+tripuls((PSDplotSpace-1)/(4*cutOffFrequency)));  

%squarer output
Outputsq = FilterOutput.^2;

%avgPSD for squaring
avgPSDsq = avgest(Outputsq,2^6);

figure
subplot(2,1,1)
plot (PSDplotSpace, abs(tPSDsq))
ylim([0 0.5])
title('Theoretical PSD of output from squarer');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot(2,1,2)
plot(linspace(0,1,N/2^6), avgPSDsq)
%ylim([0, 0.5])
title('Averaged PSD of output from squarer');
xlabel('\theta');
ylabel('Ry(\theta)');

%histogram for squarer
figure
histogram(Outputsq,100);
title('Histogram of the output from squarer');
xlabel('Amplitude');
ylabel('Occurrences');
%% theoretical and estimated PSD of half wave rectifier
%theoretical output
tPSDhwr = 1/4*(rectpuls(PSDplotSpace/(2*cutOffFrequency))+rectpuls((PSDplotSpace-1)/(2*cutOffFrequency)))+1/(4*pi)*(tripuls(PSDplotSpace/(4*cutOffFrequency))+tripuls((PSDplotSpace-1)/(4*cutOffFrequency)));

%half wave rectifier output
Outputhwr = FilterOutput;
Outputhwr(Outputhwr < 0) = 0;

%avgPSD for hwr for plotting
avgPSDhwr = avgest(Outputhwr,2^6);

figure
subplot(2,1,1)
plot(PSDplotSpace, abs(tPSDhwr))
ylim([0, 0.5])
title('Theoretical PSD of output from half wave rectifier');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot(2,1,2)
plot(linspace(0,1,N/2^6), avgPSDhwr)
ylim([0, 0.5])
title('Averaged PSD of output from half wave rectifier');
xlabel('\theta');
ylabel('Ry(\theta)');

%histogram for hwr
figure
histogram(Outputhwr,100);
title('Histogram of the output from half wave rectifier');
xlabel('Amplitude');
ylabel('Occurrences');
%% theoretical and estimated PSD of amplitude modulation
%theoretical output
f0 = 0.25;
rectshift = fftshift(tidealPSD);
tPSDam = 1/4*(circshift(rectshift,N*f0)+circshift(rectshift,-N*f0));


%amplitude modulation output
Outputam = (1 + FilterOutput).*(cos(2*pi*f0*linspace(0,N,N)))';

%avgPSD for AM for plotting
avgOutputam = avgest(Outputam, 2^6);

figure
subplot(2,1,1)
plot(PSDplotSpace, tPSDam)
ylim([0, 0.5])
title('Theoretical PSD of output from amplitude modulation');
xlabel('\theta');
ylabel('Ry(\theta)');
subplot(2,1,2)
plot(linspace(0,1,N/2^6), avgOutputam)
ylim([0, 0.5])
title('Averaged PSD of output from amplitude modulation');
xlabel('\theta');
ylabel('Ry(\theta)');

%histogram for AM
figure
histogram(Outputam,100);
title('Histogram of the output from amplitude modulation');
xlabel('Amplitude');
ylabel('Occurrences');