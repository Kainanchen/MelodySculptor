clear
close all

[source,fs] = audioread('../../Audio/SheIsMySin.wav');      % Read source
x = source(:,1);                                % Select a channel

L = size(x,1);                                  % Audio length
blockLen = fs * 5    % first 5 second window
block = x(1 : blockLen);

N = 1024;                                       % N of FFT
N2 = ceil((N+1)/2);                             % Half N
H = floor(N/2);                                 % Hop size
FH = 2^(nextpow2(0.2*fs/H)-1);                  % Frame Hop size

% F = ceil((L-N)/H);                              % No. of frames

F = ceil((blockLen-N)/H);                              % No. of frames

% X = spectrogram(x,hann(N),N-H,N,fs);            % STFT

X = spectrogram(block,hann(N),N-H,N,fs);            % STFT
X = X(round(N/40):round(N/4),:);
M = abs(X);                                     % Spectrum Magnitude
P = unwrap(angle(X));                           % Spectrum Phase

% test for negative values in v
if min(min(v)) < 0
error('matrix entries can not be negative');
return
end
if min(sum(v,2)) == 0
error('not all entries in a row can be zero');
return
end

[n,m]=size(M);
