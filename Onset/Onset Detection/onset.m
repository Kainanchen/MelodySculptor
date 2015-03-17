clear
close all
source=wavread('sheismysin');
x = source(:,1)'; 
N = length(x);

%% Compute STFT

R = 512;                
Nfft = 1024;          
X = stft(x, R, Nfft); 
    
%% Energy
m=length(X(1,:));
dse=zeros(1,m);
marke=zeros(1,m);
marked=zeros(1,m);
for i=2:m
for k=1:Nfft
    dse(i)=dse(i)+abs(X(k,i))-abs(X(k,i-1));  %energy determine function
end
delta=0.1*sum(abs(X(:,i)));  %assumed threshold
if dse(i)>delta
    marke(i)=1;
end
end
%dsem=mean(dse);
%dsemd=median(dse);
%for c=1:m;
%    if dse(c)>dsem
%        marke(c)=1;
%    end
%   if dse(c)>dsemd
%        marked(c)=1;
%    end%
%end

%% Phase
dsp=zeros(1,m);
markp=zeros(1,m);
markpd=zeros(1,m);
deltap=0.1*pi;   % random picked threshold

for i=1:m-1
    uX(:,i)=abs(X(:,i))/(sum(abs(X(:,i)))); % weighting matrix
if i<=3
    dsp(i)=0;
else
    for k=1:Nfft
        pf(k)=arg(X(k,i))-2*arg(X(k,i-1))+arg(X(k,i-2));  % phase deviation
        dsp(i)=dsp(i)+abs(princarg(pf(k)))*uX(k,i); % dsp is the weighted mean
    end
    if dsp(i)>deltap
        markp(i)=1;
    end
end
display(i);
end
%dspm=mean(dsp);
%dspmd=median(dsp);
%for d=1:m
%    if dsp(d)>dspm;
%        markp(d)=1;
%    end
%    if dsp(d)>dspmd;
%        markpd(d)=1;
%   end
%end
        
    
