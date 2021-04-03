%{
Study 2 for the course TSDT14 at Linkoping University
Improving estimates on ACF and PSD from study 1 that uses bartletts
estimate.
Author: Hoang Tran, hoatr725 (2020)
%}

clear;
clc;
close all;

%% Reconstruction of signals and estimates from study 1
%Noise construction and misc.
N = 2^16;
x = randn(N,1); %noise
ACFplotSpace = -N/2+1:1:N/2;

%low grade filt
B = ones(3,1)/3; %filter kernel
lowGradeOutput = filter(B, 1, x); %filtered with B
lgACF = bartlett(lowGradeOutput);
lgPSD = fft(lgACF);
PSDplotSpace = linspace(0,1,length(lgPSD));

tlgPSD = ((1/3 * (1+2*cos(2*pi*PSDplotSpace))).^2)';
tlgACF = circshift(fftshift(ifft(tlgPSD)),-1);

%ideal LP filt
cutOffFrequency = 0.1;
[IFB,IFA] = butter(30, cutOffFrequency*2);
idealOutput = filter(IFB, IFA, x);
idealACF = bartlett (idealOutput);
idealPSD = fft(idealACF);

tidealPSD = zeros(N,1);
tidealPSD(1:floor(cutOffFrequency*N)) = 1;
tidealPSD(ceil((1-cutOffFrequency)*N+1):N) = 1;
tidealACF = circshift(fftshift(ifft(tidealPSD)),-1);

%% Averaged PSDs

h_int = 2^8;
m_int = 2^6;
l_int = 2^4;
lgout_avgh = avgest(lowGradeOutput, h_int);
lgout_avgm = avgest(lowGradeOutput, m_int);
lgout_avgl = avgest(lowGradeOutput, l_int);
hgout_avgh = avgest(idealOutput, h_int);
hgout_avgm = avgest(idealOutput, m_int);
hgout_avgl = avgest(idealOutput, l_int);


figure
avgestplotarraylg = {abs(tlgPSD),lgout_avgh,lgout_avgm,lgout_avgl};
avgestplotarraylgtitle = {"Theoretical PSD of low grade LP filter"
                            ['Averaging over ', num2str(h_int), ' periodograms']
                            ['Averaging over ', num2str(m_int), ' periodograms']
                            ['Averaging over ', num2str(l_int), ' periodograms']};
avgestplotspacearray= {PSDplotSpace
                        linspace(0,1,N/h_int)
                        linspace(0,1,N/m_int)
                        linspace(0,1,N/l_int)};
for (i = 1:4)
    subplot(2,2,i)
    plot(avgestplotspacearray{i}, avgestplotarraylg{i});
    title(avgestplotarraylgtitle{i}); 
    xlabel('\theta');
    ylabel('Ry(\theta)');
end

figure
avgestplotarrayhg = {abs(tidealPSD),hgout_avgh,hgout_avgm,hgout_avgl};
avgestplotarrayhgtitle = {"Theoretical PSD of high grade LP filter"
                            ['Averaging over ', num2str(h_int), ' periodograms']
                            ['Averaging over ', num2str(m_int), ' periodograms']
                            ['Averaging over ', num2str(l_int), ' periodograms']};
for (i = 1:4)
    subplot(2,2,i)
    plot(avgestplotspacearray{i}, avgestplotarrayhg{i});
    title(avgestplotarrayhgtitle{i}); 
    xlabel('\theta');
    ylabel('Ry(\theta)');
end

%% windowed ACF and PSD for low grade filter
windowedlgACFrect = windowing(lgACF, @rectwin);
windowedlgACFtriang = windowing(lgACF, @triang);
windowedlgACFhamming = windowing(lgACF, @hamming);
windowedlgPSDrect = abs(fft(windowedlgACFrect));
windowedlgPSDtriang = abs(fft(windowedlgACFtriang));
windowedlgPSDhamming = abs(fft(windowedlgACFhamming));

figure
windowplotarraylg = {tlgACF
                    abs(tlgPSD)
                    windowedlgACFrect
                    windowedlgPSDrect
                    windowedlgACFtriang
                    windowedlgPSDtriang
                    windowedlgACFhamming
                    windowedlgPSDhamming};
windowplotarraylgtitle = ["Theoretical ACF of low grade LP filter"
                          "Theoretical PSD of low grade LP filter"
                          "Estimated ACF with rectangular window"
                          "Estimated PSD with rectangular window"
                          "Estimated ACF with triangular window"
                          "Estimated PSD with triangular window"
                          "Estimated ACF with hamming window"
                          "Estimated PSD with hamming window"];
for (i = 1:8)
    subplot(4,2,i)
    if (rem(i,2))
        tmp = ACFplotSpace;
        plot(tmp, windowplotarraylg{i});
        xlabel('k');
        ylabel('ry[k]');
    else
        tmp = PSDplotSpace;
        plot(tmp, windowplotarraylg{i});
        xlabel('\theta');
        ylabel('Ry(\theta)');
    end
    title(windowplotarraylgtitle{i}); 
end

figure
stemplotarraylg = [tlgACF,windowedlgACFrect,windowedlgACFtriang,windowedlgACFhamming];
stemplotarraylgtitle = [" theoretical ACF of low grade LP filter"
                            " estimated ACF with rectangular window"
                            " estimated ACF with triangular window"
                            " estimated ACF with hamming window"];
for (i = 1:4)
    subplot(4,1,i)
    stem(ACFplotSpace, stemplotarraylg(:,i));
    xlim([-20,20]);
    title(strcat('Stem plot of', stemplotarraylgtitle(i))); 
    xlabel('k');
    ylabel('ry[k]');
end

%% windowed ACF and PSD for high grade filter
windowedidealACFrect = windowing(idealACF, @rectwin);
windowedidealACFtriang = windowing(idealACF, @triang);
windowedidealACFhamming = windowing(idealACF, @hamming);

windowedidealPSDrect = abs(fft(windowedidealACFrect));
windowedidealPSDtriang = abs(fft(windowedidealACFtriang));
windowedidealPSDhamming = abs(fft(windowedidealACFhamming));

figure
windowplotarrayhg = {tidealACF
                    abs(tidealPSD)
                    windowedidealACFrect
                    windowedidealPSDrect
                    windowedidealACFtriang
                    windowedidealPSDtriang
                    windowedidealACFhamming
                    windowedidealPSDhamming};
windowplotarrayhgtitle = ["Theoretical ACF of high grade LP filter"
                          "Theoretical PSD of high grade LP filter"
                          "Estimated ACF with rectangular window"
                          "Estimated PSD with rectangular window"
                          "Estimated ACF with triangular window"
                          "Estimated PSD with triangular window"
                          "Estimated ACF with hamming window"
                          "Estimated PSD with hamming window"];
for (i = 1:8)
    subplot(4,2,i)
    if (rem(i,2))
        tmp = ACFplotSpace;
        plot(tmp, windowplotarrayhg{i});
        xlabel('k');
        ylabel('ry[k]');
    else
        tmp = PSDplotSpace;
        plot(tmp, windowplotarrayhg{i});
        xlabel('\theta');
        ylabel('Ry(\theta)');
    end
    title(windowplotarrayhgtitle{i}); 
end

figure
stemplotarrayideal = {tidealACF,windowedidealACFrect,windowedidealACFtriang,windowedidealACFhamming};
stemplotarrayidealtitle = [" theoretical ACF of high grade LP filter"
                            " estimated ACF with rectangular window"
                            " estimated ACF with triangular window"
                            " estimated ACF with hamming window"];
for (i = 1:4)
    subplot(4,1,i)
    stem(ACFplotSpace, stemplotarrayideal{i});
    xlim([-20,20]);
    title(strcat('Stem plot of', stemplotarrayidealtitle(i))); 
    xlabel('k');
    ylabel('ry[k]');
end
