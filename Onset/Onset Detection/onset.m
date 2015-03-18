clear
close all
source=audioread('sheismysin.wav');
 %source=audioread('darkchestofwonders.wav');
x = source(:,1)'; 
%y = source(:,2)'; second track
N = length(x);

%% Compute STFT

R = 256;                
Nfft = 512;  %2048 was too large for my memory        
X = stft(x, R, Nfft); 
    
%% Energy
m=length(X(1,:));
dse=zeros(1,m);
marke=zeros(1,m);
marked=zeros(1,m);
for i=2:m
%for k=1:Nfft
%    dse(i)=dse(i)+abs(X(k,i))-abs(X(k,i-1));  %energy determine function
%end
dse(i)=sum(abs(X(:,i)-abs(X(:,i-1))));
delta=0.3*sum(abs(X(:,i-1)));  %assumed threshold
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
deltap=0.5*pi;   % personally picked threshold, in the paper it stated that if there is an appearance of onset, the phase difference will be significant

for i=1:m-1
    uX(:,i)=abs(X(:,i))/(sum(abs(X(:,i)))); % weighting matrix
if i<=3
    dsp(i)=0;
else
%    for k=1:Nfft
%        pf(k)=unwrap(unwrap(angle(X(k,i))-angle(X(k,i-1)))-unwrap(angle(X(k,i-1))-angle(X(k,i-2))));  % phase deviation
%        dsp(i)=dsp(i)+abs(princarg(pf(k)))*uX(k,i); % dsp is the weighted mean
%    end
PhaseArray=unwrap([angle(X(:,i-2)),angle(X(:,i-1)),angle(X(:,i))]);
    pf=PhaseArray(:,3)-2*PhaseArray(:,2)+PhaseArray(:,1);
    dsp(i)=sum(abs(princarg(pf'*uX(:,i))));
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

%% Complex domain
disc=zeros(m);
markc=zeros(m);
H=2*ceil(m/100); % H introduces delay
for i=1:m-1
    uX(:,i)=abs(X(:,i))/(sum(abs(X(:,i))));
    if i<=3
        dsc(i)=0;
    else
%       for k=1:Nfft
%           ep(k)=princarg(2*arg(X(k,i-1))-arg(X(k,i-2)));
%            amp(k)=abs(X(k,i-1));
%            estsig=amp(k)*exp(j*ep(k));
%            disc(i)=disc(i)+sqrt((real(estsig)-real(X(k,i)))^2+(imag(estsig)-imag(X(k,i)))^2); % detection function for complex distance
%        end
        PhaseArray=unwrap([angle(X(:,i-1)),angle(X(:,i))]);
        ep=princarg(2*PhaseArray(:,2)-PhaseArray(:,1));
        amp=abs(X(:,i-1));
        estsig=amp'*exp(j*ep);
        disc(i)=sum(sqrt(real((estsig-X(:,i))).^2+(imag(estsig-X(:,i))).^2)); %This is a universal distance, should use weighted distance?
    end
        
deltac=5; % random constant?
lamda=1;
if i<H % assuming the first non-zero index is smaller than half width of median filter
    markc(i)=0;
else
    eta=median(disc(i-H/2:i+H/2));
    if disc(i-H/2)>=deltac+lamda*eta
        markc(i-H/2)=1;
    else
        markc(i-H/2)=0;
    end
end
display(i);
end;


        
        
    
