function windowedACF = windowing(ACF, windowversion)
%WINDOWING creates a windowed version of the ACF
%   Output: windowedACF     The resulting windowed ACF
%
%   Input:  ACF             The ACF to be windowed
%           windowversion   The window to be used e.g.
%                           @rectwin, @triang, @hamming etc.

N = max(size(ACF));
wind = zeros(N,1);
wind (ceil((N-65)/2):floor((N+65)/2)) = window(windowversion,65);

windowedACF = ACF.*wind;
end

