clear

U10=[1:1:51]
U10=10

M=[];
V=[];
H=[];
SS=[];
SURF=[];
US=[];
Temp10=18;
Temp0=20;
Relative_Humidity10=75;

%fetch=4.7;
fetch=5.49;
%fetch=100d3;

for i=1:length(U10)
    Uwind10=U10(i);
    
    if Uwind10<11
     ustar=sqrt(1.2d-3*Uwind10^2);
    else
    ustar=sqrt((0.49+0.065*Uwind10)*1d-3*Uwind10^2);
    end
    
fetch=100d3;


%%

Abstemp0=Temp0+273.15; %
Abstemp10=Temp10+273.15; %absolute temperature in Kelvin

g=9.81; %gravity
ka=0.4; %VonKarman Constant
Mw=18.0160d-3; %molecular weight of water kg/mol
Ms=58.443d-3; %molecular weight of NaCl kg/mol
Ma=28.9644d-3; %molecular weight of dry air kg/mol
Rg=8.314472; %gas constant J/mol/K (CODATA 99)

Salinity_w=34; %Salinity of ocean in psu
Salinity_p=34; %Salinity of particle in psu
tau_surf=(75.63 - 0.144*Temp0 + 0.221*Salinity_w)*1d-3; %Kraus and Businger 1994
rhof_init=(999.842594+6.793952d-2*Temp0-9.095290d-3*Temp0^2+1.001685d-4*Temp0^3-1.120083d-6*Temp0^4+6.536332d-9*Temp0^5); %density of fresh water kg/m^3
rhow_init=rhof_init +Salinity_w*(.824493-4.0899d-3*Temp0 +7.6438d-5*Temp0^2-8.2467d-7*Temp0^3+5.3875d-9*Temp0^4)+Salinity_w^(3/2)*(-5.72466d-3+1.0227d-4*Temp0 -1.6546d-6*Temp0^2) +4.8314d-4*Salinity_w^2; %density of salt water at 1 atm Gill p.599
rhop_init=rhof_init +Salinity_p*(.824493-4.0899d-3*Temp0 +7.6438d-5*Temp0^2-8.2467d-7*Temp0^3+5.3875d-9*Temp0^4)+Salinity_p^(3/2)*(-5.72466d-3+1.0227d-4*Temp0 -1.6546d-6*Temp0^2) +4.8314d-4*Salinity_p^2; %density of salt water at 1 atm Gill p.599
%ms=4/3*pi*R^3*rhop_init*Salinity_p/1000;
p0=101325; %pressure in Pa 

%at the surface
qs0=.62197*611.2/p0*exp((17.67*(Abstemp0-273.15))/(Abstemp0-29.66)); %Saturation humidity at surface temp and pressure from Geernaert p. 113 (Stull 1991)
rs0=qs0/(1-qs0);
r0=98/100*rs0; %Large and Pond, 1982 surface humidity over sea water (98% relative humidity at the surface)
q0=r0/(1+r0);
Tv0=Abstemp0*(1+0.6078*q0); %virtual temperature (temperature for moist air) at surface
rhoa0=1.2929*273.15/Tv0; %moist air density at the surface
nua0=1.33d-5 +(0.0084*(Tv0-273.15))*10^-5; %kinematic viscosity as function of temperature  
mua0=nua0*rhoa0; %viscosity of air
pot_T0=Abstemp0;
cpd0=1.9327d-10*Abstemp0^4-7.9999d-7*Abstemp0^3+1.1407d-3*Abstemp0^2-4.4890d-1*Abstemp0+1.0575d+03; %specific heat of dry air in J/kg/K at constant pressure
cpm0=cpd0*(1+0.84*q0); %specific heat of moist air from Makin p. 110 in Geernaert
kaT0=0.02411*(1+3.309d-3*(Abstemp0-273.15) -1.441d-6*(Abstemp0-273.15)^2); %bulk thermal conductivity from Andreas (1995) W*m/K
TDiff0=kaT0/rhoa0/cpm0; %Thermal diffusivity
Dw0=2.11d-5*(Abstemp0/273.15)^1.94*(1013.25/(p0/100)); %Diffusivity of water vapor in m^2/s Andreas(1995) (4.4) pressure=1 atm

%at 10m
p10=p0-rhoa0*g*10; %first guess pressure at 10m
pot_T10=Abstemp10; %first guess at pot_T10
for i=1:10
    qs10=.62197*611.2/p10*exp((17.67*(pot_T10-273.15))/(pot_T10-29.66)); %saturation humidity at 10 m temp and pressure from Geernaert p. 113 (Stull 1991)
    rs10=qs10/(1-qs10);
    r10=Relative_Humidity10/100*rs10; %calculating 10 m r from relative humidity and saturation mixing ratio
    q10=r10/(1+r10);%calculating 10 m q from relative humidity and saturation humidity
    
    Tv10=pot_T10*(1+0.6078*q10);%virtual temperature (temperature for moist air) at 10 m
    rhoa10=1.2929*273.15/Tv10; %moist air density at 10m
    p10=p0-0.5*(rhoa10+rhoa0)*g*10; %pressure at 10m
    ke10=2/7*((1-q10+q10/.62197)/(1-q10+8*q10/7/.62197));%Gill
    pot_T10=Abstemp10*(p0/p10)^ke10; %ams glossary and gill ----alternative by P&D as function of z only
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of AIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% SURFACE WAVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_m=sqrt(rhow_init*g/tau_surf); %Bond waves
    X=fetch; %fetch in m
    k_0=g/Uwind10/Uwind10; %plants lower wind-wave wavenumber limit

    X=X*k_0; %normalized fetch
    X0=2.2d4;
    Omega_c=0.84*(tanh((X/X0)^0.4))^-0.75; %wave age parameter
    theta=(-pi/2:pi/100:pi/2); %wave-wind angles -pi/2 to pi/2  index 51 is 0
    dtheta=gradient(theta); %
    wlength=logspace(log10(2*pi/1500),4.5,1500);
    wnumber=2*pi./wlength; %wave number
    wspeed=sqrt(g./wnumber.*(1+(wnumber/k_m).^2)); %wave speed in m/s 
    omega=wnumber.*wspeed;
    frequency=omega./2/pi;
    deltaf=abs(gradient(frequency));
    deltak=abs(gradient(wnumber));
    [Sk, B_sat, k_peak, p_spr]=myunified_spectrum_s(wnumber, theta, Uwind10,ustar,Omega_c,k_m);
    p_spr=2*p_spr;
    
    Sf=B_sat.*wnumber.^(-3).*deltak./deltaf;
    Spk=B_sat.*wnumber.^(-3).*deltak;
    wheight=sqrt(abs(Spk)); %amplitude of wave
    
    peak=find(Spk==max(Spk));
    omegapeak=omega(peak);
    frequencypeak=frequency(peak);
    periodpeak=1/frequencypeak;
    wlengthpeak=wlength(peak);
    p_spr=min(p_spr, p_spr(peak));
    
    significantwheight=sqrt(sum(abs(Spk)))*4;
%%
     
q=q10; %humidity (specific) 
qs=qs10; %Saturation humidity at air temp
f=Relative_Humidity10/100; %Relative humidity/100
AbsTemp=Abstemp10;
Droptemp=Abstemp10;
Tv=Tv10;
nua=1.33d-5 +(0.0084*(Tv-273.15))*10^-5;
Pressure=p10;
Ur=Uwind10;
rhow=rhof_init; %density of pure water
rhop=rhop_init; %density of pure water
rhoa=rhoa10;
S=Salinity_p;

%% Charateristic
%Estimate of inertial time scale and terminal velocity;

%  esat=611.2*exp((17.67*(AbsTemp-273.15))/(AbsTemp-29.66)); %from Fairall
%  rhoma=1.2929*273.13/Tv; %density of moist air  %generates lower evaporation temperature 
%  rhovapor=q*rhoma; % water vapor density
%  rhoad=1.2929*273.13/AbsTemp; %density of dry air 
% 
 Lv=(25.00895-.02274*(AbsTemp-273.15))*10^5; %Latent heat of evaporation from Andreas (1995) (2.4) Fleagle and Businger (1980) in J/kg
%  Dw=2.11d-5*(AbsTemp/273.15)^1.94*(1013.25/(Pressure/100)); %Diffusivity of water vapor in m^2/s Andreas(1995) (4.4) pressure=1 atm
%  deltaw=8d-8; %in m empirical constant from Andreas (1995) p. 857
%  alphac=0.036; %empirical constant from Andreas (1995) p. 857
%  alphac=1; %empirical constant from Lewis and Schwartz p59
% 



%Dprimew=Dw./(R./(R + deltaw) + Dw./R./alphac*(2*pi*Mw/Rg/AbsTemp)^(1/2)); %Diffusivity of water vapor for noncontinuum effects 
% 
 cps=4217.4 -3.720283*(Droptemp-273.15)+0.1412855*(Droptemp-273.15)^2-2.654387d-3*(Droptemp-273.15)^3 +2.093236d-5*(Droptemp-273.15)^4;%specific heat of water from Gill appendix 3 in J/kg/K for fresh water
 cps=cps +S*(-7.6444+0.107276*(Droptemp-273.15)-1.3839d-3*(Droptemp-273.15)^2)+S^(3/2)*(0.17709-4.0772d-3*(Droptemp-273.15)+5.3539d-5*(Droptemp-273.15)^2); %specific heat of sea water droplet in J/kg/K
% 

%Dprimew=Dw./(R./(R + deltaw) + Dw./R./alphac*(2*pi*Mw/Rg/AbsTemp)^(1/2)); %Diffusivity of water vapor for noncontinuum effects 
 cpd=1.9327d-10*AbsTemp^4-7.9999d-7*AbsTemp^3+1.1407d-3*AbsTemp^2-4.4890d-1*AbsTemp+1.0575d+03; %specific heat of dry air in J/kg/K at constant pressure
%  deltaT=2.16d-7; % in m empirical constant from Andreas (1995)
%  alphat=0.7; %empirical constant from Andreas (1995)
%  ka=0.02411*(1+3.309d-3*(AbsTemp-273.15) -1.441d-6*(AbsTemp-273.15)^2); %bulk thermal conductivity from Andreas (1995) W*m/K
%  kprimea=ka./(R./(R + deltaT) + ka./R./alphat/rhoad/cpd*(2*pi*Ma/Rg/AbsTemp)^(1/2)); 
 
 
%  delta=Droptemp/AbsTemp -1; %from Andreas (1990) p. 482
%  nu=2; %number of ions NaCl dissociates into
%  phis=0.91154614+1.7317496707d-4*S+4.7616058412d-6*S^2-9.2541509027d-9*S^3+7.3475024678d-12*S^4; %osmotic coefficient good for S>30
%  sigmas=235.8*(1-Droptemp/647.096)^1.256*(1-0.625*(1-Droptemp/647.096))*1d-3+0.221*S*1d-3; %surface tension incuding temperature and salinity effects
%  esat_drop=611.2*exp((17.67*(Droptemp-273.15))/(Droptemp-29.66));
% % 
%  nua=1.33d-5 +(0.0084*(Tv-273.15))*10^-5; %kinematic viscosity of moist air
%  fv=1+0.25*sqrt(2*R*Ur/nua); %ventilation coefficient for vapor diffusion Meirink p.102
%  fh=fv; %ventilation coefficient for heat diffusion
%  ms=4/3*pi*R.^3*rhop*S/1000;
%  y=(2*Mw*sigmas/Rg/AbsTemp/(1+delta)/rhow./R - nu*phis*ms*(Mw/Ms)./((4*pi*R.^3*rhop/3)-ms));
% 
%  dr=fv./R.*Dprimew*Mw*esat/rhop/Rg/AbsTemp.*(f-1/(1+delta)*exp(Lv*Mw/Rg/AbsTemp*(delta/(1+delta)))*exp(y)); %modified Andreas (1990)/Edson and Fairall (1994)
%  Qr=Mw*esat_drop/Rg/Droptemp*exp(y); %vapor density at droplet surface Andreas (1990) p. 483
%  Qa=rhovapor; % water vapor density
%  dT=-3./R.^2/rhop/cps.*(fh.*kprimea.*(AbsTemp-Droptemp) +fv.*Lv.*Dprimew.*(Qa-Qr)); %modified Andreas (1990) /Edson and Fairall 1994
%  dr=1/rhow./R.*fv.*Dprimew.*(Qa-Qr); %To be consitent with above; not assuming esat_drop=esat*exp(Lv*Mw/Rg/AbsTemp) Pruppacher 13.25 (1997)
 
 

%%%%%%%%%%%%%% RADIUS dist%%%%%%%%%%%%
R=logspace(-6,-3,100); R=R';

  
 
[Vt_p, Cd_p, RE, CU, Tau_I, Tau_stokes]=vp4(R, rhop, rhoa, nua);

%% Distribution of Fairall

r_80=0.518*(R*1d6).^0.976;
for ttt=1:length(R)
if (r_80(ttt))>0.9 && (r_80(ttt))<15
        dFdr_F(ttt)=75.7672/2*3.8d-6*Uwind10^3.4*(r_80(ttt))^-0.024*10^(4.405-2.646*log(r_80(ttt))/log(10)^(1)-3.156*log(r_80(ttt))^2/log(10)^(2)+8.902*log(r_80(ttt))^3/log(10)^(3)-4.482*log(r_80(ttt))^4/log(10)^(4));
    elseif (r_80(ttt))>=15 && (r_80(ttt))<=37.5
        dFdr_F(ttt)=75.7672/2*3.8d-6*Uwind10^3.4*(r_80(ttt))^-0.024*1.02d4*(r_80(ttt))^(-1);
    elseif (r_80(ttt))>=37.5 && (r_80(ttt))<=100
        dFdr_F(ttt)=75.7672/2*3.8d-6*Uwind10^3.4*(r_80(ttt))^-0.024*6.95d6*(r_80(ttt))^(-2.8);   
    elseif (r_80(ttt))>100
        dFdr_F(ttt)=75.7672/2*3.8d-6*Uwind10^3.4*(r_80(ttt))^-0.024*1.75d17*(r_80(ttt))^(-8);
    else   dFdr_F(ttt)=NaN; 
    end    
end


  U14=Uwind10+ustar/ka*log((14)/(10));
       A98(1)=10.0^(0.0676*U14 + 2.43);
       A98(2)=10.0^(0.959*sqrt(U14) - 1.476);
% %	  FIND CONSTANTS FOR EXTRAPOLATING TO LARGE DROPLETS.
       DFDR80_10 = A98(1)*exp(-3.1*((log(10/2.1))^2))+ A98(2)*exp(-3.3*((log(10/9.2))^2));
% %	  THIS IS THE FUNCTION FOR THIS WIND SPEED AT R80 = 10.
       C1 =10.0*DFDR80_10;
       C2 =(37.5^1.8)*C1;
       C3 =(100.0^5.2)*C2;
%        
%     %%for 15 m/s
ratio=(A98(1)*exp(-3.1*((log(r_80(max(find((r_80)<10)))/2.1))^2))+ A98(2)*exp(-3.3*((log(r_80(max(find((r_80)<10)))/9.2))^2)))/(3.5*(C1*(r_80(min(find((r_80)>10))))^-1)*0.506*(r_80(min(find((r_80)>10))))^-0.024);
for ttt=1:length(R)
    if (r_80(ttt))>0.9 && (r_80(ttt))<10
        dFdr_A(ttt)= A98(1)*exp(-3.1*((log(r_80(ttt)/2.1))^2))+ A98(2)*exp(-3.3*((log(r_80(ttt)/9.2))^2))/ratio;
     elseif (r_80(ttt))>10 && (r_80(ttt))<37.5
         dFdr_A(ttt)=3.5*(C1*(r_80(ttt))^-1)*0.506*(r_80(ttt))^-0.024;
     elseif (r_80(ttt))>=37.5 && (r_80(ttt))<=100
         dFdr_A(ttt)=3.5*(C2*(r_80(ttt))^-2.8)*0.506*(r_80(ttt))^-0.024;
     elseif (r_80(ttt))>100
         dFdr_A(ttt)=3.5*(C3*(r_80(ttt))^-8)*0.506*(r_80(ttt))^-0.024;
     else   dFdr_A(ttt)=NaN; 
     end  
end



for ttt=1:length(R)
if (r_80(ttt))>0.5 && (r_80(ttt))<0.9
indx=sum(r_80<0.9)+1;
a=dFdr_F(indx)*r_80(indx)^3;
dFdr_F(ttt)=a*r_80(ttt)^-3; 
end
end

%[Vt_p, Cd_p, RE, CU, Tau_I, Tau_stokes]=vp4(R, rhop, rhoa, nua);


dFdr=dFdr_F'*1d6; %in meters!
Vd_p=Vt_p./(1-exp(-Vt_p./1d-3/Uwind10)); 
%Vd_p=Vt_p;%includes turbulence (smith et al 93)
dCdr=dFdr./Vd_p;
dR=gradient(R); 

r_bar=nansum(dCdr.*R.*dR)/nansum(dCdr.*dR);
dr=sqrt(nansum(dCdr.*(R-r_bar).^2.*dR)/nansum(dCdr.*dR));
phi_volume=nansum(4/3*pi*dCdr.*R.^3.*dR);
phi_mass=nansum(rhop/rhoa*4/3*pi*dCdr.*R.^3.*dR);
phi_sensibleheat=nansum(cps/cpd*rhop/rhoa*4/3*pi*dCdr.*R.^3.*dR);
phi_totalheat=nansum((cps*Droptemp+Lv)/(cpd*Droptemp)*rhop/rhoa*4/3*pi*dCdr.*R.^3.*dR);

volume_flux=nansum(4/3*pi*dFdr.*R.^3.*dR);
Prod_rate=nansum(dFdr.*dR);


M= [M phi_mass];
V=[V nansum(4/3*pi*dFdr.*R.^3.*dR)];
SS=[SS Prod_rate];
H=[H phi_totalheat];
SURF=[SURF nansum(4*pi*dFdr.*R.^2.*dR)];
US=[US ustar];


end

 figure
 loglog(R*1d6, dFdr/1d6)
 hold on
 axis([1 1000 1d-4 1d6]) 
 
 figure
 loglog(R*1d6, 4/3*pi*dFdr.*R.^3)
 hold on
 axis([1 1000 1d-4 1d6])
 
  figure
 loglog(R*1d6, dCdr/1d6) % plot in microns
 hold on
 axis([1 1000 1d-1 1d6]) 
 

%% Choix du rayon representatif

V_o=Vd_p;
U_o=U10;
eta=V_o./U_o;
R_o=R;
x_o=significantwheight;
m_o=rhow_init*4/3*pi*R_o.^3;
%[Vt_o, Cd_p_o, RE_o, CU_o, Tau_Io, Tau_stokes_o]=vp4(Ro, rhop, rhoa, nua);
t_o=x_o./V_o;
dC_odr=dCdr;

%dans le paper - gradeur characteristiques
gamma=7/5;
Cv=cpm0/gamma;
Cps=cps;
rho_sw=rhow_init; %density of sea water
rho_o=rhoa10; %density of air
T_o=Temp10; %temperature of air
mu_o=nua*rho_o; %viscosity of air
Lv=Lv; %Latent heat of evaporation
ka_o=0.02411*(1+3.309d-3*(Abstemp10-273.15) -1.441d-6*(Abstemp10-273.15)^2);



yo=15.14; a=5.6e-4; b=1.338; %RH= 75%
%yo=17.59; a=6.7e-5; b=1.364; %RH= 95%
T_ev=yo+a*exp(-b*log10(R));%RH= 75%
R_eq=0.46*R; alpha=0.46; %RH= 75%
%R_eq=0.713*R; %RH= 95%
u_eq=U10;

Tau_Io=Tau_I;
Tau_To=10.^(-4.343088735+1.6066803348*log10(R.*1d6)+0.1297515062*(log10(R.*1d6)).^2-0.0146703809*(log10(R.*1d6)).^3);
 
Tau_R_75RH=10.^(-1.4555786169+1.5485214917*log10(R.*1d6)+0.1260301567*(log10(R.*1d6)).^2-0.0133380536*(log10(R.*1d6)).^3);
Tau_R_95RH=10.^(-0.6069852307+1.6697599561*log10(R.*1d6)+0.0999532367*(log10(R.*1d6)).^2-0.0135869637*(log10(R.*1d6)).^3);
a=((f-0.95)/(0.75-0.95))^1/1; 
Tau_Ro=Tau_R_75RH*a+(1-a)*Tau_R_95RH;

%%
T=logspace(-5,5,100);%time
DT=gradient(T);

%NEW
for jjt=1:length(T) %loop in time
    m0=4/3*pi*R.^3*rhow_init;
    
    %density
    zeta=T(jjt)./Tau_Ro; Phi=(1-alpha)*(3*alpha^2*exp(-zeta)+3*alpha*(1-alpha)*exp(-2*zeta)+(1-alpha)^2*exp(-3*zeta));
    MASS(jjt)=nansum(4/3*pi*R.^3*rhow.*(1-alpha^3-Phi).*dCdr.*dR); %That's density really!!!
    HUMIDITY(jjt)=(1-q)/rhoa*MASS(jjt);
    
    
    %Velocity
    zeta=T(jjt)./Tau_Io;
    U(jjt)=-1/rhoa*U10*alpha^2*(1+rhow/rhow_init)*nansum(m0.*(1-exp(-zeta./alpha^2)).*dCdr.*dR);
    
    %Sensible Heat
    zeta=T(jjt)./Tau_To;
    TS(jjt)=cps/(rhoa*Cv)*nansum((m0-6*alpha^2*(1-alpha)*4/3*pi*R.^3*rhow).*(1-exp(-zeta)).*(Temp0-T_ev).*dCdr.*dR);
    
     %Latent Heat
    TL(jjt)= -Lv/(rhoa*Cv)*nansum(4/3*pi*R.^3*rhow.*(1-alpha^3-Phi).*dCdr.*dR);
end

QQ=(q+HUMIDITY)./(1-(q+HUMIDITY))./rs10; %relative humidity after spray

figure
semilogx(T,QQ*100)

figure
semilogx(T,U+U10)
figure
semilogx(T,TS+Temp10)
figure
semilogx(T,TL+Temp10)

s=[QQ'*100 U'+U10 TS'+Temp10' TL'+Temp10];
save D:\people\Mieussens\Notes_on_Kinetic_Theory\s.dat s -ascii

MASS_LIMIT=(1-alpha^3)*nansum(4/3*pi*R.^3*rhow.*dCdr.*dR);
U_LIMIT=-1/rhoa*U10*alpha^2*(1+rhop/rhow_init)*nansum(4/3*pi*R.^3*rhow.*dCdr.*dR);

%%
%OLD
for jj=1:length(R);
rayon=R(jj);
mass_goutte=4/3*pi*rayon.^3*rhop;
surface_goutte=4*pi*rayon.^2;
SENS(jj,:)=cps*mass_goutte/Tau_To(jj)*(Temp0-T_ev(jj))*exp(-T/Tau_To(jj)); %Sensible for each drop as a finction of time
LAT(jj,:)=Lv*surface_goutte*rhow/Tau_Ro(jj)*(R(jj)-R_eq(jj))*exp(-T/Tau_Ro(jj)); %Latent heat for each drop as a finction of time
DENSITY(jj,:)=surface_goutte*rhow/Tau_Ro(jj)*(R(jj)-R_eq(jj))*exp(-T/Tau_Ro(jj)); %mass for each drop as a finction of time
MOMENTUM(jj,:)=mass_goutte/Tau_Io(jj)*(0-u_eq)*exp(-T/Tau_Io(jj));
end

 s=[SENS(1,:)' SENS(34,:)' SENS(67,:)' SENS(end,:)' LAT(1,:)' LAT(34,:)' LAT(67,:)' LAT(end,:)' DENSITY(1,:)' DENSITY(34,:)' DENSITY(67,:)' DENSITY(end,:)' MOMENTUM(1,:)' MOMENTUM(34,:)' MOMENTUM(67,:)' MOMENTUM(end,:)' ];
save D:\people\Mieussens\Notes_on_Kinetic_Theory\s.dat s -ascii

%integral in radius accounting for the concentration
for jjt=1:length(T)
    SENS_Total(jjt)=nansum(SENS(:,jjt).*dCdr.*dR); %J/s=Watt (rho*Cv*dT/dt) 
    LAT_Total(jjt)=nansum(LAT(:,jjt).*dCdr.*dR); %J/s=Watt (rho*Cv*dT/dt)
    DENSITY_Total(jjt)=nansum(DENSITY(:,jjt).*dCdr.*dR); %density - it's pi_m !!!!
    MOMENTUM_Total(jjt)=nansum(MOMENTUM(:,jjt).*dCdr.*dR);
end

Total_Temp_Induced_Sensible=nansum(SENS_Total.*DT)/(rhoa*Cv); %K
Total_Temp_Induced_Latent=nansum(LAT_Total.*DT)/(rhoa*Cv);
Total_density_Induced=nansum(DENSITY_Total.*DT);%
Total_humidity_Induced=Total_density_Induced*(1-q)/rhoa;%
Total_momentum_Induced=nansum(MOMENTUM_Total.*DT);%
Total_speed_Induced=nansum(MOMENTUM_Total.*DT)/rhoa;%


TS=Temp10+cumsum(SENS_Total.*DT)/(rhoa*Cv); %K
TL=Temp10+cumsum(LAT_Total.*DT)/(rhoa*Cv); %temperature after spray
semilogx(T,TS)
semilogx(T,TL)

Q=cumsum(DENSITY_Total.*DT)*(1-q)/rhoa;
QQ=(q+Q)./(1-(q+Q))./rs10; %relative humidity after spray
semilogx(T,QQ*100)

SP=U10+cumsum(MOMENTUM_Total.*DT)/rhoa;  %wind speed after spray
semilogx(T,SP)


s=[QQ'*100 SP' TS' TL'];
save D:\people\Mieussens\Notes_on_Kinetic_Theory\s.dat s -ascii
%%



%rapports de constantes:
Cps/Cv;
rho_sw/rho_o;
Lv/(Cv*T_o);

%particules
%droplet_radius=Ro;  
Tau_mo=rho_sw/(3*rhow)*Tau_Ro;
Tau_go=V_o/g;


phi_o=1/rho_o*m_o.*dC_odr;

Ma=U_o/sqrt(gamma*Pressure/rho_o); % Mach number
Re=rho_o*U_o*x_o/mu_o; %Reynolds number (airflow at scale x_0)
Pe=cpm0*x_o*rho_o*U_o/ka_o;
Fr=U_o/sqrt(x_o*g);

%rapports de constantes dans l'equation de WB


alpha_I=t_o./Tau_Io;
alpha_R=t_o./Tau_Ro;
alpha_m=t_o./Tau_mo;
alpha_T=t_o./Tau_To;
alpha_g=t_o./Tau_go;


%Pour William Boltzman
figure
loglog(R, alpha_I,'k'); hold on
loglog(R, alpha_m,'g')
loglog(R, alpha_T,'b')
loglog(R, alpha_R,'r')
legend('alpha_I','alpha_m','alpha_T', 'alpha_R')

figure
loglog(R, Tau_Io,'k'); hold on
loglog(R, Tau_mo,'g')
loglog(R, Tau_To,'b')
loglog(R, Tau_Ro,'r')
loglog(R, t_o,'o')

%rapports de constantes dans l'equation de la masse:
A=1./eta; % transport
theta_m=alpha_m.*phi_o; %A2=alpha_mphi; %perte de masse

%rapports de constantes dans l'equation de la vitesse:
AA=1./eta; % transport
BB=1./(gamma*Ma.^2.*eta); %pression
CC=1./(Re.*eta); %diffusion
DD=1./(Fr.^2.*eta); %gravitee 
theta_D=alpha_I.*phi_o;  %trainee
theta_m=alpha_m.*phi_o;  %perte de momentum parceque perte de masse


%Figure pour l'equation de la vitesse
figure
loglog(R, theta_D,'o'); hold on
loglog(R, theta_m,'sq')
loglog(R, CC,'<')
loglog(R, DD,'>')
loglog(R,BB,'*')
legend('Drag', 'mass', 'Reynolds (diffusion)', 'Froude (gravity)', 'Pressure')



%rapports de constantes dans l'equation de la chaleur:
theta_T1=Cps/Cv*phi_o; %transport (mixture) phi_o is Xi in paper
AAA=1./eta; % transport
CCC=1./eta.*gamma./Pe; %diffusion (heat from from Q)
DDD=1./eta.*gamma*(gamma-1)*Ma.^2./Re; %friction air-air
%EEE=1./eta.*(gamma-1)*Ma.^2./Re; %friction air-air
EEE=1./eta.*(gamma-1); %pressure

theta_T2=Cps/Cv*theta_m; %perte de quantite de chaleur sensible parceque perte de masse
theta_T3=theta_m; %perte de quantite de chaleur sensible parceque perte de masse;  theta_m=Alpha_m*Xi in paper
theta_T4=Lv/(Cv*T_o)*theta_m; %perte de chaleur latente 
theta_T5=gamma*(gamma-1)*Ma.^2.*theta_m; %perte de chaleur  parceque perte d'energie cinetique
theta_T6=Cps/Cv*alpha_T.*phi_o; %variation de temperature (sensible)
theta_T7=gamma*(gamma-1)*Ma.^2.*theta_D; %travail de la force de trainee

%Figure pour l'equation de la chaleur
%close all
figure 
loglog(R,theta_T2,'ok'); hold
loglog(R,theta_T3,'sq')
loglog(R,theta_T4,'*')
loglog(R,theta_T5,'.')
loglog(R,theta_T6,'x')
loglog(R,theta_T7,'v')
loglog(R,CCC,'<')
legend( 'theta_T2 transport','theta_T3 chaleur sensible masse','theta_T4 chaleur latent masse', 'theta_T5 cinetique','theta_T6 sensible heat','theta_T7 trainee', 'CCC Peclet (diffusion)', 'DDD Reynolds (friction air)')




figure
loglog(R,Ma,'o')
hold on
loglog(R,1./un_sur_FR_2.*Ma,'*')
legend('Mach','Froude carre sur Mach')





