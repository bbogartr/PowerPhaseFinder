fs=250;
T=1;
L=fs*T;
t=linspace(0,1,fs);
f=11;
A=5;
w=2*pi*f;
% make random phases

% make random signals  with f signal embedded in random noise
n=3;
Sigs=zeros(n,L);
for j=1:n
    a=rand(1)*2*pi;
    Sigs(j,:)=(rand(1)*1)+sin((w*t)+a)+randn(1,L);
%     sb=(rand(1)*1)+sin((w*t)+b)+randn(1,fs);
%     sc=(rand(1)*1)+sin((w*t)+c)+randn(1,fs);
end
% Sigs=[sa;sb;sc];

% % look at signals
% figure;
% plot(sa)
% hold on
% plot(sb,'-g')
% plot(sc,'-k')

% spectrogram(sa,[],[],[],250) 


phz=zeros(n,1);
pow=zeros(n,1);
figure

pltrow=ceil(n/2);
for j=1:n
%fourier
[Y]=fft(Sigs(j,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
% plot(P1(2:end))
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of S(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

alpha_peak=find(P1(2:end)==max(P1(2:end)));
% use floor to round to nearest whole radian
% s_peak_phz=floor(angle([Y(alpha_peak+1)]));
phz(j,1)=angle([Y(alpha_peak+1)]);
pow(j,1)=P1(alpha_peak+1);
theta=repmat(phz(j,1),10);
rho=linspace(0,pow(j,1),10);

subplot(pltrow,pltrow,j)
polarplot(theta,rho,'-r','LineWidth',2)
rlim([0 1])
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(['Signal ', num2str(j)])
% soundsc(s)
end


% dist=zeros(1,((n^2)-n)/2);
dist=zeros(n);
% x=0;
for j=1:n
    for k=1:n
%         x=x+1;
        dist(j,k)=phz(j)-phz(k);
    end
end
dist=dist./(2*pi);
dist=tril(dist,-1);

dist_tri=[];
for j=1:n
    dist_tri=[dist_tri,dist(j+1:end,j)'];
end


figure;
%show similarity matrix
subplot(1,2,1)
imagesc(squareform(dist_tri))
title("difference in phases between signals (green is 0)")
%lower threshold to abs<.1 (chosen as a "low" sounding number) and find "chord" 
thresh=abs(dist_tri)<.1;

subplot(1,2,2)
imagesc(squareform(thresh))
title("'chord' phase diff below .1")




%% Make notes


mfs=14400;
rpow=round(pow*10)/10
rpow=[.5,.7,.3]
%fs for each electrode from notes file
f1=notes(max(find(notes(:,2)<rpow(1) & notes(:,1)==1))+1,3);
f2=notes(max(find(notes(:,2)<rpow(2) & notes(:,1)==2))+1,3);
f3=notes(max(find(notes(:,2)<rpow(3) & notes(:,1)==3))+1,3);

%three seconds for each note
T=1;
L=mfs*T;
t=linspace(0,T,L);

w1=2*pi*f1;
w2=2*pi*f2;
w3=2*pi*f3;

% make random phases

% make tone signals
n=3;
Sings=zeros(n,L);
Sings(1,:)=(thresh(1)|thresh(2))*sin((w1*t));
Sings(2,:)=(thresh(1)|thresh(3))*sin((w2*t));
Sings(3,:)=(thresh(2)|thresh(3))*sin((w3*t));

%only play each note if it is synched with another note
Song=Sings(1,:)+Sings(2,:)+Sings(3,:);
player = audioplayer(Song, 14400);
play(player)