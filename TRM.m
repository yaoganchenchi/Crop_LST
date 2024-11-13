function [output] = TRM(SWin, LWin, Ta, q, albedo, ra, rs, emis, rhoa, Ps, G)

globalconstant= getConstants();
mask = ~isnan(SWin) & ~isnan(LWin) & ~isnan(Ta) & ~isnan(q)...
& ~isnan(albedo) & ~isnan(ra) & ~isnan(rs) & ~isnan(emis) & ~isnan(rhoa) & ~isnan(Ps) & ~isnan(G);
SWin(~mask)=nan;
LWin(~mask)=nan;
Ta(~mask)=nan;
q(~mask)=nan;
albedo(~mask)=nan;
ra(~mask)=nan;
rs(~mask)=nan;
emis(~mask)=nan;
rhoa(~mask)=nan;
Ps(~mask)=nan;
G(~mask)=nan;

%% deal with ground heat flux and its derivatives with respect to albedo, ra and rs
Rn_star = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ta.^4;

%% important intermediate variables
delta   = desdT(Ta);

gamma = (globalconstant.cp.*Ps)./(globalconstant.epsilon.*globalconstant.Lv); %(globalconstant.cp.*Ps)./(0.622.*globalconstant.Lv);

[qsat,~]    = qs(Ta, Ps);

lambda0 = 1./(4*globalconstant.sb.*emis.*(Ta.^3));

ro = rhoa.*globalconstant.cp.*lambda0;

ftrm = ro./ra.*(1+delta./gamma.*(ra./(ra+rs)));

AA = rhoa.*globalconstant.Lv.*(qsat-q);

% intermedia derivatives
dftrm_dRA = -ro./(ra.^2).*(1+delta./gamma.* (ra./(ra+rs)).^2);

dftrm_dRS = -delta./gamma.*ro./(ra+rs).^2;

dRn_star_dTA = -4*globalconstant.sb.*emis.*(Ta.^3);

dQ_sat_Ta_dTA = 0.622./ (Ps) .*delta; % check it

dR0_dTA = -rhoa.*globalconstant.cp.*(3/4)./(globalconstant.sb.*emis.*(Ta.^4));

ddelta_dTA =des2dT2(Ta);

dftrm_dTA = dR0_dTA./ra .* (1+delta./gamma.*ra./(ra+ra)) + ro./ra .* ddelta_dTA./gamma.*ra./(ra+rs);

dlambda0_dTA = -3/4 ./(globalconstant.sb.*emis.*(Ta.^4));


%% TRM sensitivities
output.dTs_dALBEDO = -lambda0.*SWin./(1+ftrm);
output.dTs_dRA = lambda0.*AA./(ra+rs).^2./(1+ftrm)  - dftrm_dRA.*lambda0.*(Rn_star-G-AA./(ra+rs))./(1+ftrm).^2;
output.dTs_dRS = lambda0.*AA./(ra+rs).^2./(1+ftrm)  - dftrm_dRS.*lambda0.*(Rn_star-G-AA./(ra+rs))./(1+ftrm).^2;
output.dTs_dEMIS = lambda0.*( LWin-globalconstant.sb.*(Ta.^4) - ((Rn_star-G) - AA./(ra+rs))./(emis.*(1+ftrm)))./(1+ftrm); % corrected on Aug 11
output.dTs_dSWin = lambda0.*(1-albedo) ./(1+ftrm);
output.dTs_dLWin = lambda0.*emis ./(1+ftrm);
output.dTs_dQA = rhoa.*globalconstant.Lv./(ra+rs)./(1+ftrm).*lambda0;
output.dTs_dTA = lambda0.*(dRn_star_dTA-rhoa.*globalconstant.Lv./(ra+rs).*dQ_sat_Ta_dTA)./(1+ftrm) + dlambda0_dTA.*(Rn_star-G-AA./(ra+rs))./(1+ftrm) - dftrm_dTA.*lambda0 .*(Rn_star-G-AA./(ra+rs))./(1+ftrm).^2 + 1;
output.dTs_dG = -lambda0./(1+ftrm);


%% computation of Ts, H, LE, and energy balance based on linearized SEB but radiatively coupled (this is from Raupach 2001, QJRMS and Rigden/Li 2017, GRL)

Ts = Ta + lambda0.*(Rn_star - G -AA./(ra+rs))./(1+ftrm);

H = rhoa.*globalconstant.cp.*(Ts - Ta)./ra;

[Qsat_LE,~] = qs(Ts, Ps);
LE = rhoa.*globalconstant.Lv.*(Qsat_LE-q)./(ra+rs);

energy_balance = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4  - H - LE - G;

Rn = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4;

output.Ts=Ts;
output.H=H;
output.LE=LE;
output.energy_balance=energy_balance;
output.Rn=Rn;

end

function delta2= des2dT2(T)
delta2 = 611*exp(((1727*T)./100 - 9434601/2000)./(T - 717/20)).*(1727./(100*(T - 717/20)) - ...
    ((1727*T)/100 - 9434601/2000)./(T - 717/20).^2).^2 - 611*exp(((1727*T)./100 - ...
    9434601/2000)./(T - 717/20)).*(1727./(50*(T - 717/20).^2) - (2*((1727*T)./100 - 9434601/2000))./(T - 717/20).^3);

end

function delta= desdT(T1)
% T in K, delta is Pa/K
delta = 611.*exp(((1727*T1)/100 - 9434601/2000)./(T1 - 717/20)).*(1727./(100*(T1 - 717/20)) - ((1727*T1)/100 - 9434601/2000)./(T1 - 717/20).^2);
end

function [ qsat,esat ] = qs( T, P )
globalconstant= getConstants();
esat = 611*exp(17.27*(T-273.15)./(T-273.15+237.3));  % Eq. 3.9a from Dingman in Pa
qsat = globalconstant.epsilon*esat./P; % not account for vapor pressure 
end

function globalconstant= getConstants()
globalconstant.cp  = 1004.64;      % specific heat at constant pressure, J/kg/K
globalconstant.Lv  = 2.4665*10^6; % latent heat of vaparization

globalconstant.R   = 287.058;       % dry air gas constant J/kg/K
globalconstant.Rv  = 461.5;       % water vapor gas constant J/kg/K
globalconstant.g   = 9.80616;      % gravity constant m/s^2
globalconstant.k   = 0.4 ;      % von-Karman constant
globalconstant.p0  = 101325 ;      % pa standard
globalconstant.emis  = 0.95 ;             % emissivity, for veg
globalconstant.sb          = 5.670367*10^(-8); % stephan-boltzman constant, W/(m^2 K^4)
globalconstant.Avogadro = 6.02214*10^26; % avogadro's number
globalconstant.k_bol = 1.38065*10^-23; % boltzmann constant
globalconstant.MWda = 28.966; % molecular weight of dry air

globalconstant.lapse_rate = globalconstant.g/globalconstant.cp;

globalconstant.epsilon= globalconstant.R/globalconstant.Rv  ; 
globalconstant.gamma = globalconstant.cp*globalconstant.p0./(globalconstant.epsilon*globalconstant.Lv);  % psychrometric constant at p0

end













