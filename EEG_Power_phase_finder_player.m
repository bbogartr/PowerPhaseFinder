% EEG = pop_loadset('filename','br_ICApruned.set','filepath','D:\\EEGdata\\TutorialCleaningdata\\4_ICA_pruned\\');
% EEG=load('E:/EEG_summer/EEGdata/DotComp_br.easy')
load('A:/Documents/EEG_music/matlabscripts/guitarnotes.mat')

fs=500; %sampling rate (like that of EEG)
T=1;%how long the signal is in seconds
L=fs*T; %how long the signal is in samples
t=linspace(0,T,fs); %t for making signals 
f=12; %frequency to make (11 Hz, alpha)
w=2*pi*f;
% make random phases
bigsong=zeros(1,length(EEG.data));
beg=1;
% powees=zeros(n,round(length(EEG.data)/fs));
% phzees=zeros(n,round(length(EEG.data)/fs));
cc=0
for z=1:fs:(length(EEG.data)-(L*5))
    % make random signals  with f signal embedded in random noise
    cc=cc+1;
    beg=beg+L;
    n=3;
    Sigs=zeros(n,L);
    chans=[4,7,11];
    for j=1:n
        
%         a=rand(1)*2*pi;%random phase between 0 and 2pi to add below
        %sigs(j,:)=radnom amplitude 0:1+sinewave at random phase + random noise
    %     Sigs(j,:)=(rand(1))+sin((w*t)+a)+randn(1,L);
        Sigs(j,:)=EEG.data(chans(j),beg:beg+L-1);
        
    end
    %get phase and power for each electrode signal
    phz=zeros(n,1);
    pow=zeros(n,1);
    for j=1:n
    %fourier
    [Y]=fft(Sigs(j,:));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    psdx = (1/(fs*L)) * abs(P1).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    alpha_peak=find(P1(2:end)==max(P1(2:end)));
    % use floor to round to nearest whole radian
    % s_peak_phz=floor(angle([Y(alpha_peak+1)]));
    phz(j,1)=angle([Y(alpha_peak+1)]);
%     phzees(j,cc)=phz(j,1);
    pow(j,1)=P1(alpha_peak+1)/4;
%     powees(j,cc)=pow(j,1);
    end
    % dist=zeros(1,((n^2)-n)/2);
    dist=zeros(n);
    % x=0;
    for j=1:n
        for k=1:n
            dist(j,k)=phz(j)-phz(k);
        end
    end
    dist=dist./(2*pi);
    %get the bottom traingle of matrix, dist
    dist=tril(dist,-1);

    dist_tri=[];
    for j=1:n
        dist_tri=[dist_tri,dist(j+1:end,j)'];
    end

    %lower threshold to abs<.1 (chosen as a "low" sounding number) and find "chord" 
    thresh=abs(dist_tri)<.02;
    %% Make notes

    mfs=14400;
    rpow=round(pow*10)/10;
    %fs for each electrode from notes file
    %This treats electrode as a string and the power as the fret
    %less than .2 power is open, .2 power and up is 1st, 2nd ,3rd fret
    f1=notes(find(notes(:,2)<rpow(1) & notes(:,1)==1,1,'last'),3);
    f2=notes(find(notes(:,2)<rpow(2) & notes(:,1)==2,1,'last'),3);
    f3=notes(find(notes(:,2)<rpow(3) & notes(:,1)==3,1,'last'),3);

    %three seconds for each note
    T=1;
    L2=mfs*T;
    t=linspace(0,T,L2);

    % make tone signals
    %only play each note if it is synched with another note (use logical values from thresh)
    Sings=zeros(n,L2);

    w1=2*pi*f1;
    w2=2*pi*f2;
    w3=2*pi*f3;

    Sings(1,:)=(thresh(1)|thresh(2))*sin((w1*t)); % if electrodes 1&2 | 1&3 are synch = 1
    Sings(2,:)=(thresh(1)|thresh(3))*sin((w2*t));
    Sings(3,:)=(thresh(2)|thresh(3))*sin((w3*t));
    %add them together
    Song=Sings(1,:)+Sings(2,:)+Sings(3,:);
    bigsong=[bigsong,Song];
end

player = audioplayer(bigsong, 14400);
play(player)