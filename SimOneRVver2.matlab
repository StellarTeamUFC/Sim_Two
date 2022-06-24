clear all
close all

% Planetary and solar parameters
Msun=1.989E+30; % in kg
Mearth=5.972E+24; % in kg
Mjupiter=1.899E+27 % in kg

% input parameters
Porb=50; % Orbital period in days
i=15; % inclination angle of star (between 0 and 90 degrees)
omega=42; % argument of pericenter in true orbital plane
omega=pi/180*omega % in rad
tp=1; % in days

np=inputdlg('Planet Mass fraction (in Jupiter mass units):','PMfrac');
np=str2num(np{:});

ns=inputdlg('Stellar Mass fraction (in solar mass units):','SMfrac');
ns=str2num(ns{:});


e=0.9; % eccentricity
Mp=np*Mjupiter; % planet mass
Mstar=ns*Msun; % stellar mass

% Radial velocity semi-amplitude (K)
K=203.29.*((Mp.*sin(pi/180*i))./Mjupiter).*((Msun./(Mstar+Mp)).^(2/3)).*((Porb).^(-1/3)).*(1./sqrt(1-e^2)); % in m/s

% Calculate of the true anomaly; angle referred to the elliptical (true) orbit.

j=0;

for t=1:100
j=j+1;
syms Et

eqnLeft=(2*pi/Porb)*(t-tp); % mass function M(t)

eqnRight=Et-e.*cos(Et); % our goal: find Et

EtS=vpasolve(eqnLeft == eqnRight,Et);
EtS=double(EtS);
EtSj(j)=double(EtS);

% time vector in days
time(j)=tp+(Porb/(2*pi)).*(EtS-e.*sin(EtS)); % new equation of time

% Calculate the ni(t)
nit=acos((cos(EtS)-e)./(1-e.*cos(EtS)));
nitj(j)=nit;

% Radial velocity
RV(j)=K.*(cos(omega+nit)+e.*cos(omega));

end

RV=RV';
time=time';

plot(time,RV,'k.');
xlabel('time (days)','FontSize',12);
ylabel('Radial Velocity (m/s)','FontSize',12);

% Include: noisy time series (white noise)
stdnoise=5; % noise SD
fs=length(time)-1;
noise=ffgn(stdnoise,0.500001,1,fs+1,0)'; %H=0.5 (Gaussian noise)

figure;

RVnoise=RV+noise; % noisy RV time series
plot(time,RVnoise,'k.');
xlabel('time (days)','FontSize',12);
ylabel('Radial Velocity (m/s)','FontSize',12);

figure;
% Remove periodically points (2 each 3 points)

RVnoiseRed=downsample(RVnoise,3);
timeRed=downsample(time,3);

plot(timeRed,RVnoiseRed,'.');
xlabel('time (days)','FontSize',12);
ylabel('Radial Velocity (m/s)','FontSize',12);

figure;
EtSj=EtSj';
plot(time,EtSj,'.');
xlabel('time (days)','FontSize',12);
ylabel('E(t)','FontSize',12);

figure;
nitj=nitj';
plot(time,nitj,'.');
xlabel('time (days)','FontSize',12);
ylabel('ni(t)','FontSize',12);



