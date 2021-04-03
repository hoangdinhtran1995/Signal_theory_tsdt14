function estimatedACF = bartlett(signal)
%BARTLETT uses the bartlett estimate to give an approximate result for the
%ACF of a filter
%   output: estimatedACF is the bartlett estimated ACF
%   input:  signal is the signal to have its ACF estimated

N = size(signal);
N = N(1);

estimatedACF = zeros(N,1);

for(k = -N/2+1:N/2)

    sum = 0; %set sum initially to 0    

    for (n = 1:(N-abs(k))) %calc sum for given position k
        sum = sum+signal(n+abs(k))*signal(n);
    end
        
    estimatedACF (N/2+k) = 1/N*sum;
end
end

