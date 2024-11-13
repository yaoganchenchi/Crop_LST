function A004_prepare_TRM_inputs_yearly_diurnal_MODIS_ERA(iyear)
if ~exist('iyear','var')
    iyear = 2010;
end
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
addpath(sprintf('%scode/TRM/',root_dir));

ipt.path_in = sprintf('/projectsp/f_cc2127_1/chenchi/P001_SUHII/ERA5_land/'); %ERA5;
ipt.path_in_albedo = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/MCD43C3_extracted/';
ipt.path_in_LST = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/LST/extracted/MOD11C3/all/gap_filled/';
ipt.path_in_EMIS = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/LST/extracted/MOD11C3/all/gap_filled/';
ipt.path_in_SWin = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/MCD18C1_extracted/';
ipt.path_in_LE = '/projectsp/f_cc2127_1/chenchi/Data_Archive/LE_C61_MOD_composite/CMG_0.05deg_monthly/';
ipt.path_in_LAI = '/projectsp/f_cc2127_1/chenchi/Data_Archive/LAIFPAR_C61_MCD_composite/CMG_0.05deg_monthly/';


path_out = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_input_eramodis/yearly/diurnal/');
if ~(exist(path_out,'dir') ==7)
    mkdir(path_out);
end

tic
% constants
constants=getConstants(); % get global constants
constants.flux_scale = 86400;
constants.syear = iyear;
constants.eyear = iyear;
ipt.R = georefcells([-90 90],[-180 180],[1800,3600],'ColumnsStartFrom','north');
ipt.row=3600;
ipt.col=7200;
ipt.nmonth=12;

%% 1 & 2 get LH and SH annual
get_sshf_slhf(ipt,path_out,constants)

%% 3 Emissivity
get_EMIS(ipt,path_out,constants)
%
%% 4 & 5 2m air temperature, skin temperature,
get_2t_skt(ipt,path_out,constants)

%%  6 & 7  & 8 Shortwave in,  Longwave in, Albedo,
get_ssrd_strd(ipt,path_out,constants)
get_fal(ipt,path_out,constants)

%% 9 specific humidity
get_Qair(ipt,path_out,constants)

%% 10 & 11  air density, air pressure
get_RHOA(ipt,path_out,constants)
get_sp(ipt,path_out,constants)

%% 12 ground heat flux
get_GH(ipt,path_out,constants)

%% 13 & 14 aerodynamic resistance, surface resistance
get_RA(ipt,path_out,constants) 
get_RS(ipt,path_out,constants) 

toc
end

%% SH_LH
function get_sshf_slhf(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

%% SH
save_root = sprintf('%ssshf/',path_out);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

data_root = sprintf('%ssshf/',ipt.path_in);
for iyear=constants.syear:constants.eyear
    fprintf('ERA5 monthly: sshf, %04d\n',iyear)
    fn = sprintf('%sERA5_land_%s_%04d.mat',data_root,'sshf',iyear);
    data=loadMatData(fn);
    new_data=mean(data,3,'omitnan')/constants.flux_scale;
    output=georesize(new_data,ipt.R,2,'bilinear');
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'sshf',iyear);
    save(savename,'output','-v7.3');
    clear output
end

%% LH
save_root2 = sprintf('%sslhf/',path_out);
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end
data_root = sprintf('%s',ipt.path_in_LE);
for iyear=constants.syear:constants.eyear
    fprintf('MODIS monthly: slhf, %04d\n',iyear)
    data=nan(ipt.row,ipt.col,ipt.nmonth,'single');
    for imonth=1:ipt.nmonth
        fn = sprintf('%sMOD_LE_C61_0.05CMG_monthly_%04d%02d_filtered.mat',data_root,iyear,imonth);
        temp = loadMatData(fn);
        data(:,:,imonth)=temp.LE;
    end
    output=-mean(data,3,'omitnan');
    savename = sprintf('%s%s_%04d_yearly.mat',save_root2,'slhf',iyear);
    save(savename,'output','-v7.3');
    clear output

end

end

%% Emissivity
function get_EMIS(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

%% EMIS
save_root = sprintf('%sEMIS/',path_out); % annual
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

emis_dir = sprintf('%s',ipt.path_in_EMIS);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5 monthly: EMIS, %04d\n',iyear)
    emis = loadMatData(sprintf('%sMOD11C3.005CMG.filled.A%04d.EMIS.mat',emis_dir,iyear));
    
    % save annualy
    output = mean(emis,3,'omitnan');
    output(output>1)=1;
    output(output==0)=nan;
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'EMIS',iyear);
    save(savename,'output','-v7.3');
    clear output
end

end

%% skt temp and 2m air temperature
function get_2t_skt(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

save_root2 = sprintf('%s2t/',path_out);
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end

save_root3 = sprintf('%sskt/',path_out);
if ~(exist(save_root3,'dir') ==7)
    mkdir(save_root3);
end

skt_dir = sprintf('%s',ipt.path_in_LST);
t2_dir =sprintf('%s2t/',ipt.path_in);
emis_dir = sprintf('%s',ipt.path_in_EMIS);
d2_dir =sprintf('%s2d/',ipt.path_in);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5 2t, MODIS skt, %04d\n',iyear)
    emis = loadMatData(sprintf('%sMOD11C3.005CMG.filled.A%04d.LST.mat',emis_dir,iyear));
    skt = loadMatData(sprintf('%sMOD11C3.005CMG.filled.A%04d.LST.mat',skt_dir,iyear));
    d2 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',d2_dir,'2d',iyear));
    t2 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',t2_dir,'2t',iyear)) - 9.8*98/1000;
    t2(t2<d2)=d2(t2<d2);
    
    % 2t
    data = (mean(constants.sb.*t2.^4,3,'omitnan') ./constants.sb) .^0.25 ;
    output=georesize(data,ipt.R,2,'bilinear');
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root2,'2t',iyear);
    save(savename,'output','-v7.3');
    clear output
    
    %skt
    step1 = mean(emis.*skt.^4,3,'omitnan');
    step2 = mean(emis,3,'omitnan');
    output = (step1./step2) .^0.25 ;
    savename = sprintf('%s%s_%04d_yearly.mat',save_root3,'skt',iyear);
    save(savename,'output','-v7.3');
    clear output

end
end

%% Shortwave in, longwave in
function get_ssrd_strd(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
save_root = sprintf('%sssrd/',path_out);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

save_root2 = sprintf('%sstrd/',path_out);
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end

ssrd_dir = sprintf('%s',ipt.path_in_SWin);
strd_dir = sprintf('%sstrd/',ipt.path_in);

for iyear=constants.syear:constants.eyear
    fprintf('MODIS ssrd, ERA strd, monthly %04d\n',iyear)
    
    ssrd = nan(ipt.row,ipt.col,ipt.nmonth,'single');
    for imonth=1:ipt.nmonth
        ssrd(:,:,imonth) = loadMatData(sprintf('%sMCD18C1.A%04d%02d.mat',ssrd_dir,iyear,imonth));
    end
    
    strd = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',strd_dir,'strd',iyear));
    
    % ssrd
    output=mean(ssrd,3,'omitnan');
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'ssrd',iyear);
    save(savename,'output','-v7.3');
    clear output
    
    %strd
    data=mean(strd,3,'omitnan')/constants.flux_scale;
    output=georesize(data,ipt.R,2,'bilinear');
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root2,'strd',iyear);
    save(savename,'output','-v7.3');
    clear output
end


end

%% ALBEDO
function get_fal(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
save_root = sprintf('%sfal/',path_out);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

fal_dir = sprintf('%s/',ipt.path_in_albedo);
for iyear=constants.syear:constants.eyear
    fprintf('MODIS monthly: fal, %04d\n',iyear)
    
    fal = nan(ipt.row,ipt.col,ipt.nmonth,'single');
    for imonth=1:ipt.nmonth
        fal(:,:,imonth) = loadMatData(sprintf('%sMCD43C3.A%04d%02d.mat',fal_dir,iyear,imonth));
    end
    
    % fal
    output=mean(fal,3,'omitnan');
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'fal',iyear);
    save(savename,'output','-v7.3');
    clear output

end

end

%% specific humidity
function get_Qair(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
addpath(sprintf('%scode/TRM/',root_dir));

save_root = sprintf('%sQair/',path_out);  % yearly
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

save_root2 = sprintf('%sQair/',ipt.path_in); % monthly
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end

d2_dir = sprintf('%s2d/',ipt.path_in); % 2m dew point temp
sp_dir = sprintf('%ssp/',ipt.path_in);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5 monthly: Qair, %04d\n',iyear)
    d2 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',d2_dir,'2d',iyear));
    sp = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',sp_dir,'sp',iyear));
    
    % Qair, monthly
    [output,~]=qs(d2,sp);
    
    savename = sprintf('%sERA5_land_%s_%04d.mat',save_root2,'Qair',iyear);
    save(savename,'output','-v7.3');
    clear output
    
    % Qair, yearly
    [temp,~]=qs(d2,sp);
    data=mean(temp,3,'omitnan');
    output=georesize(data,ipt.R,2,'bilinear');
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'Qair',iyear);
    save(savename,'output','-v7.3');
    clear output

end

end

%% air density RHOA
function get_RHOA(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

save_root = sprintf('%sRHOA/',path_out); % yearly
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

save_root2 = sprintf('%sRHOA/',ipt.path_in); % monthly
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end

d2_dir = sprintf('%s2d/',ipt.path_in); % 2m dew point temp
sp_dir = sprintf('%ssp/',ipt.path_in);
t2_dir =sprintf('%s2t/',ipt.path_in);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5 monthly: RHOA, %04d\n',iyear)
    sp = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',sp_dir,'sp',iyear));
    t2 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',t2_dir,'2t',iyear))- 9.8*98/1000;
    d2 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',d2_dir,'2d',iyear));
    [~,ea]=qs(d2,sp);
    t2(t2<d2)=d2(t2<d2);
    
    % RHOA, monthly
    output = (sp-ea)./ (constants.R.*t2) + ea./ (constants.Rv.*t2);
    savename = sprintf('%sERA5_land_%s_%04d.mat',save_root2,'RHOA',iyear);
    save(savename,'output','-v7.3');
    clear output
    
    % yearly
    data = mean((sp-ea)./ (constants.R.*t2) + ea./ (constants.Rv.*t2),3,'omitnan');
    output=georesize(data,ipt.R,2,'bilinear');
    
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'RHOA',iyear);
    save(savename,'output','-v7.3');
    clear output
end

end

%% air pressure
function get_sp(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

save_root = sprintf('%ssp/',path_out); % yearly
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

sp_dir = sprintf('%ssp/',ipt.path_in);
for iyear=constants.syear:constants.eyear
    fprintf('ERA5 monthly: sp, %04d\n',iyear)
    sp = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',sp_dir,'sp',iyear));
    % yearly
    data = mean(sp,3,'omitnan');
    output=georesize(data,ipt.R,2,'bilinear');
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'sp',iyear);
    save(savename,'output','-v7.3');
    clear output
end
end

%% ground heat flux
function get_GH(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

save_root = sprintf('%sGH/',path_out); % yearly
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

save_root2 = sprintf('%sGH_monthly/',path_out); % monthly
if ~(exist(save_root2,'dir') ==7)
    mkdir(save_root2);
end

ssrd_dir = sprintf('%s',ipt.path_in_SWin);
strd_dir = sprintf('%sstrd/',ipt.path_in);
slhf_dir = sprintf('%s',ipt.path_in_LE);
sshf_dir = sprintf('%ssshf/',ipt.path_in);
fal_dir = sprintf('%s',ipt.path_in_albedo);
EMIS_dir = sprintf('%s',ipt.path_in_EMIS);
skt_dir = sprintf('%s',ipt.path_in_LST);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5+MODIS monthly: GH, %04d\n',iyear)
    strd0 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',strd_dir,'strd',iyear));
    sshf0 = loadMatData(sprintf('%sERA5_land_%s_%04d.mat',sshf_dir,'sshf',iyear));
    EMIS = loadMatData(sprintf('%sMOD11C3.005CMG.filled.A%04d.EMIS.mat',EMIS_dir,iyear));
    skt = loadMatData(sprintf('%sMOD11C3.005CMG.filled.A%04d.LST.mat',skt_dir,iyear));
    
    ssrd = nan(ipt.row,ipt.col,ipt.nmonth,'single');
    slhf = nan(ipt.row,ipt.col,ipt.nmonth,'single');
    fal = nan(ipt.row,ipt.col,ipt.nmonth,'single');
    strd= nan(ipt.row,ipt.col,ipt.nmonth,'single');
    sshf= nan(ipt.row,ipt.col,ipt.nmonth,'single');
    for imonth=1:ipt.nmonth
        ssrd(:,:,imonth) = loadMatData(sprintf('%sMCD18C1.A%04d%02d.mat',ssrd_dir,iyear,imonth));
        fal(:,:,imonth) = loadMatData(sprintf('%sMCD43C3.A%04d%02d.mat',fal_dir,iyear,imonth));
        
        fn = sprintf('%sMOD_LE_C61_0.05CMG_monthly_%04d%02d_filtered.mat',slhf_dir,iyear,imonth);
        temp = loadMatData(fn);
        slhf(:,:,imonth)=-temp.LE;
        
        
        strd(:,:,imonth) = georesize(strd0(:,:,imonth),ipt.R,2,'bilinear');
        sshf(:,:,imonth) = georesize(sshf0(:,:,imonth),ipt.R,2,'bilinear');
    end
    
    
    % monthly
    GH = ssrd.*(1- fal) + EMIS.*strd./constants.flux_scale...
        - EMIS.*constants.sb.*skt .^4 +sshf ./ constants.flux_scale + slhf;
    output = GH;
    savename = sprintf('%s%s_%04d.mat',save_root2,'GH',iyear);
    save(savename,'output','-v7.3');
    clear output
    
    % yearly
    output = mean(GH,3,'omitnan');
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'GH',iyear);
    save(savename,'output','-v7.3');
    clear output
end

end


%% aerodynamic resistance RA - average meteorology first 
function get_RA(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
save_root = sprintf('%sRA/',path_out);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end

sshf_dir = sprintf('%ssshf/',path_out);
skt_dir = sprintf('%sskt/',path_out);
t2_dir = sprintf('%s2t/',path_out);
RHOA_dir = sprintf('%sRHOA/',path_out);
for iyear=constants.syear:constants.eyear
    fprintf('ERA5 + MODIS monthly: RA, %04d\n',iyear)
    sshf = loadMatData(sprintf('%s%s_%04d_yearly.mat',sshf_dir,'sshf',iyear));
    skt = loadMatData(sprintf('%s%s_%04d_yearly.mat',skt_dir,'skt',iyear));
    t2 = loadMatData(sprintf('%s%s_%04d_yearly.mat',t2_dir,'2t',iyear));
    RHOA = loadMatData(sprintf('%s%s_%04d_yearly.mat',RHOA_dir,'RHOA',iyear));
    
    % yearly
    RA = RHOA.*constants.cp.*(skt-t2)./(-sshf);
    mask = RA<0 | RA==-Inf | RA==Inf;
    RA(mask)=nan;
    
    output = RA;
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'RA',iyear);
    save(savename,'output','-v7.3');
    clear output  
end

end

%% surface resistance RS  -   average meteorology first 
function get_RS(ipt,path_out,constants)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

save_root = sprintf('%sRS/',path_out);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end


skt_dir = sprintf('%sskt/',path_out);
RHOA_dir = sprintf('%sRHOA/',path_out);
Qair_dir = sprintf('%sQair/',path_out);

slhf_dir = sprintf('%sslhf/',path_out);
sp_dir = sprintf('%ssp/',path_out);
RA_dir = sprintf('%sRA/',path_out);

for iyear=constants.syear:constants.eyear
    fprintf('ERA5 + MODIS monthly: RS, %04d\n',iyear)
    slhf = loadMatData(sprintf('%s%s_%04d_yearly.mat',slhf_dir,'slhf',iyear));
    sp = loadMatData(sprintf('%s%s_%04d_yearly.mat',sp_dir,'sp',iyear));
    
    skt = loadMatData(sprintf('%s%s_%04d_yearly.mat',skt_dir,'skt',iyear));
    Qair = loadMatData(sprintf('%s%s_%04d_yearly.mat',Qair_dir,'Qair',iyear));
    RHOA = loadMatData(sprintf('%s%s_%04d_yearly.mat',RHOA_dir,'RHOA',iyear));
    RA = loadMatData(sprintf('%s%s_%04d_yearly.mat',RA_dir,'RA',iyear));
    
    [Qs_sat,~]=qs(skt,sp);
    
    Lv = nan(size(sp));
    Lv(:,:) = constants.Lv;  
    
    RS = RHOA./(-slhf./Lv).*(Qs_sat-Qair)-RA;
    mask = RS<0  | RS==-Inf | RS==Inf;
    RS(mask)=nan;
    
    %yearly
    output = RS;
    savename = sprintf('%s%s_%04d_yearly.mat',save_root,'RS',iyear);
    save(savename,'output','-v7.3');
    clear output
end

end