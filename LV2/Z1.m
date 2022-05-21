clear all; close all; clc;
s = tf('s');
w = 0.01:0.001: 110;
Kc = 190;
G = Kc*0.1/((0.2*s+1)*(0.02*s+1))

[mag,phase] = bode(G, w);
mag = squeeze(mag);
phase = squeeze(phase);

figure(1)
semilogx(w, 20*log10(mag))
text(22.33, interp1(w, 20*log10(mag), 22.33), sprintf('\\leftarrow  dB at %.2f rad/sec',22.33), 'HorizontalAlignment','left', 'VerticalAlignment','middle')
grid on
%%
wq = interp1(w, 20*log10(mag), 22.33)
%%
wq = interp1(20*log10(mag), w,  -10*log10(1/alpha))
%%
zc = 0.517;
pc = 1.77;
Gc = (s/zc + 1) / ( s/pc + 1);
Gtotal = Gc * G;
bode(G)
hold on; grid on;
bode(Gtotal)

%%
Gp = 1/(s^2 + 1);
step(Gp/(1+Gp))

%%
s = tf('s');
G = 144000/(s*(s+36)*(s+100));
[mag,phase,wout] = bode(G);
mag = squeeze(mag);
phase = squeeze(phase);
wq = interp1(20*log10(mag), wout, -3.76);               % Find Desired Frequency
figure
semilogx(wout, 20*log10(mag), wq, -3.76, 'r+')
grid
text(wq, -3.76, sprintf('\\leftarrow -3.76 dB at %.2f rad/sec',wq), 'HorizontalAlignment','left', 'VerticalAlignment','middle')


%%

clear all; close all; clc;
s = tf('s');
w = 0.01:0.001: 110;
PMspec = 50; % zeljena fazna margina
Gp = 2/((s+1)*(s+2)*(s+3)); % proces

% 1.korak => rucno odredimo broj integratora i pojacanje kc na osnovu e(t) i dobijemo
% novi sistem koji posmatramo G
kc = 2.5;  
G = Gp * kc / s;

% 2.korak => očitamo fazu nekompenziranog sistema i wx
[mag,phase] = bode(G, w); % ne vrati u dB (samo plota u dB), phase je u stepenima
mag = squeeze(mag);
phase = squeeze(phase);
wx = interp1(20*log10(mag), w, 0) % tamo gdje je 20log10(mag) = 0
phaseGwx = interp1(w, phase, wx); % tamo gdje je w = wx
PMuncom = 180 + phaseGwx % faza nekompenziranog sistema

% 3. korak  => odredimo FImax i alpha
FImax = PMspec - PMuncom + 10
radijan = pi / 180;
alpha = (1 - sin(FImax * radijan)) / (1 + sin(FImax*radijan))

% 4. korak => trazimo wMax
wMax = interp1(20*log10(mag), w, -10*log10(1/alpha))

% 5.korak => odraditi pc i zc jer je wMax = sqrt(pc*zc) i alpha = zc/pc

zc = wMax * sqrt(alpha)
pc = zc / alpha


Gc = (s/zc + 1) / ( s/pc + 1);
Gtotal = Gc * G;
bode(Gtotal, w)

%%
s=tf('s');
a=18;
G=(3.75*s+a*0.3)/(s*(10*s+2)*(s+1)*(s+1.5+a));
t=0:0.01:10;
nagib=1; %jedinicna rampa
rampa=nagib*t; 


%syms w; 
%G_p_f=(3.75i*w+a*0.3)/((1i*w)*(10i*w+2)*(1i*w+1)*(1i*w+1.5+a));
%w180=double(solve(phase(G_p_f)==pi, w));
 %T_krit=abs(2*pi/w180) 
 %K_krit=1/abs(freqresp(G, w180)) 
 %Kp=K_krit; Ki=0;
 
 
 %Kp=0.35*K_krit; 
 %Ki=Kp/(1.25*T_krit); 
 %Kd=Kp*(0.2*T_krit);
 %G_reg= Kp + Ki*(1/s);
 %G_reg=Kp+ Ki*(1/s)+ Kd*s;
 %C=G_reg;

%pidTuner(G)
Kp=6.235;
Kd=20.96;
C=Kp+Kd*s;

%za info na step
Gk=C*G/(1+C*G)
stepinfo(Gk)

%pidTuner(G)

bode(Gk)
margin(Gk)
step(Gk)

b=1;
tau=a+b+10;
%Guk=G*exp(-tau*s);
lamda=1;
n=3;
C_imc=1/(G*(lamda*s+1)^n)