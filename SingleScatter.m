clear all; close all
%%
% This code is implemented to investigate the single scattering solutions
% and the "colours" of the light reaching the observer from an horizontal
% path. The atmosphere may contain molecules and aerosol.
% Definition of parameters that don't change in the code
% g         is Henyey-Greenstein g parameter
% bt        is Brightness Temperature (BT)

g=0.85;
bt=5800;    %BT of the Sun

% Definition of parameters that change in the code
% zena      is the sun zenith angle
% scata     is the scattering angle
% maratio   is the ratio of molecules to aerosol scattering at 0.56 microns
% alpha     is the Angstrom coefficient

% Different cases are valuated expressed by job='A#', where # is the number
% of case analyzed

job='A0'; zena=40; scata=130; maratio=0; alpha=1;
% job='A1'; zena=40; scata=130; maratio=1; alpha=1; 
% job='A2'; zena=40; scata=130; maratio=1; alpha=3; 
% job='A3'; zena=40; scata=130; maratio=10; alpha=1;
% job='A4'; zena=70; scata=160; maratio=10; alpha=1; 
% job='A5'; zena=40; scata= 50; maratio=10; alpha=1;
%%
% Choice of wavelength to compute VISibility and normalise particle 
% scattering coefficient

wd=0.01;                        % Step in wavelength vector
iwa=27;                         % wa(iwa)= 0.52 microns

%wd=0.02; 
%iwa=6;

wa=[0.3:wd:1.0]; 
n=size(wa,2);                   % Wavelength in micron
mu=cos(scata*pi/180);           % Conversion from degrees to radians
mus=cos(zena*pi/180);           % Conversion from degrees to radians

for i=1:n
% How to compute Rayleigh optical depth tauf(i) for a vertical path through 
% the atmosphere with airmass=1
% Frohlich and Shaw, AO, 1980 vol 19 n. 11, 1773-1775
barg=3.916+0.074*wa(i)+0.05/wa(i);
tauf(i)=8.38e-3/wa(i)^barg;
end

Hv=7.7;                         % Vertical scale height: same for molecular
                                % scattering coefficient ksm and aerosol 
                                % scattering coefficient ksa

ksm=tauf./Hv;                   % Rayleigh molecular scattering coefficient    
ksa=ksm;                        % Rayleigh aerosol scattering coefficient

ksm(i)=1.1e-3 / wa(i)^4.09;     % Origin unknown

% Normalize aerosol scattering coefficient using maratio
k_aernorm=ksm(iwa)*maratio/wa(iwa)^(-alpha);

for i=1:n
ksa(i)=k_aernorm * wa(i)^(-alpha);
end

d=25;                           % Distance from observer to black wall in m
kst=ksm + ksa;                  % Total scattering coefficient
tau_v=Hv*kst;                   % Total vertical optical depth for a path
                                % of length d
tau_h=d*kst;                    % Total horizontal optical depth for a path
                                % of length d
% Connection with visibility
vis=3.912/(ksm(iwa)+ksa(iwa));
B=1.e-6;                        % Assumes a value of A=0.20 m2
% Definition of parameters for Planck function
c0=2.997*10.^8;                 % m/s speed of light in vaccum
h=6.625*10.^-34;                % J.s Planck constant 
k_=1.38*10.^-23;                % T/K Boltzmann constant
nn=1;                           % Refravtive index of the medium
for i = 1:n
planck =(2*pi.*h.*(c0.^2))./((nn.^2).*(wa(i).^5).*(exp((h.*c0)./(nn.*k_.*bt.*wa(i).*i))-1));
Ts(i)=exp(-tau_v(i)/mus);
E0(i)=planck * B * mus;         % mus to compute irradiance on a horiz.path
E0S(i)=E0(i)*Ts(i);
end

Pmol=0.750*(1 + mu*mu);             % Scattering direction 
Paer=(1-g^2)/(1+g^2-2*g*mu)^(3/2);  % Paer is assumed independent of wavelength
Pt=(Pmol*ksm + Paer*ksa)./(ksm+ksa);
Th=exp(-tau_h);

% Definition of Radiance L
L=E0S .* Pt .* (1-Th);          % E0S is irradiance [Wm-2]; Pt [1/sr];

% Normalization of L, E0 and E0S

LN=L/max(L);
E0N=E0/max(E0);
E0SN=E0S/max(E0S);

% Definition of axis and elements on the graphics
% The function num2str(A) converts a numeric array into a character array 
% that represents the numbers

st0=num2str(d);
st1=num2str(maratio);
st2=num2str(vis,'%7.1f');
st3=num2str(alpha);

% Concatenate two or more character vectors or two or more cells array
str1=strcat('d=',st0,',','M/A=',st1,',',' \alpha =',st3,' Vis=',st2);

st0=num2str(zena); 
st1=num2str(scata);
str0=strcat('-',job,' sz=',st0,' sa=',st1);

% Figures
figure(1);
subplot(3,2,1);
hold on;
set(gca,'XMinorTick','on','YMinorTick','on', 'box','on');
hold on
plot(wa,ksm,'-b', 'linewidth',2);
plot(wa,ksa,'-r', 'linewidth',2);
H = legend('ksm','ksa');
xlabel('Wavelength [micron]');
ylabel('Scattering coeff.[1/m]');
title(str0);

subplot(3,2,2);
hold on;
set(gca,'XMinorTick','on','YMinorTick','on', 'box','on');
hold on;
plot(wa,tau_v,'-b', 'linewidth',2);
plot(wa,tau_h,'-r', 'linewidth',2);
H = legend('tau_v','tau_h');
xlabel('wavelength (micron)'); 
ylabel('Optical Path'); 
title(str1);

subplot(3,2,3);
set(gca,'XMinorTick','on','YMinorTick','on', 'box','on');
hold on;
plot(wa,Ts,'-b', 'linewidth',2);
plot(wa,Th,'-r', 'linewidth',2);
H = legend('Ts','Th','Location','NW');
xlabel('Wavelength [micron]');
ylabel('Transmittance');
title(str1);

subplot(3,2,4);
rad=pi/180;
ang=[0:0.5:10 11:5:180 180];
na=length(ang); 
mua=cos(ang*rad);
was=[1, floor(n/2), n];     % The function floor(A) round matrix elements 
                            % toward negative infinity 
sa=length(was);
for j=1:na
Pm(j)=0.750*(1 + mua(j)*mua(j));        % Molecular scattering
Pa(j)=(1-g^2)/(1+g^2-2*g*mua(j))^(3/2); % Aerosol scattering
for k=1:sa
is=was(k);
Ptot(k,j)=(Pm(j)*ksm(is) + Pa(j)*ksa(is))./(ksm(is)+ksa(is)); % Tot. Scatt.
end
end

% How to convert axis in semi-log
semilogy(ang,Pm,'--b','Linewidth',3);
hold on;
semilogy(ang,Pa,'--r','Linewidth',3);
semilogy(ang,Ptot(1,:),'-k','Linewidth',2);
semilogy(ang,Ptot(2,:),'-g','Linewidth',2);
semilogy(ang,Ptot(3,:),'-m','Linewidth',2);
H = legend('Pm','Pa','Ptot(0.30)','Ptot(0.64)','Ptot(1.0)','Location','SW');
xlabel('Scattering angle');
ylabel('Scattering diagram');

subplot(3,2,5);
set(gca,'XMinorTick','on','YMinorTick','on', 'box','on');
hold on;
plot(wa,E0,'-g','linewidth',2);
plot(wa,E0S,'-b', 'linewidth',2);
plot(wa,L,'-r', 'linewidth',2);
H = legend('E0','E0S','L');
xlabel('wavelength (micron)'); 
ylabel('Irradiance');
title(str1);

subplot(3,2,6);
set(gca,'XMinorTick','on','YMinorTick','on', 'box','on');
hold on;
plot(wa,E0N,'-g', 'linewidth',2);
plot(wa,E0SN, '-b','linewidth',2);
plot(wa,LN,'-r', 'linewidth',2);
H = legend('E0N','E0SN','LN','Location','SE');
xlabel('wavelength (micron)');
ylabel('Relative units');
title(str1);






