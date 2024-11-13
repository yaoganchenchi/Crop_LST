function A005_prepare_TRM_inputs_yearly_diurnal_MODIS_ERA_adj(iyear)
if ~exist('iyear','var')
    iyear = 2000;
end
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
addpath(sprintf('%scode/TRM/',root_dir));

fn_lc = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/output_c61/MCD12C1.CMG005.C61.IGBP.MODE.LC.2001.2021.mat';
fn_srcwin = sprintf('%ssearch_window_size_MCD12.mat','/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/');

lc = loadMatData(fn_lc);
srcwin = loadMatData(fn_srcwin);

% read dem
dem_fn = sprintf('%sdata/GMTED2010_DEM_CMG005.tif',root_dir);
dem = readgeoraster(dem_fn);
dem_save_root = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_input_eramodis/yearly/diurnal_nearbyall/dem/';
if ~(exist(dem_save_root,'dir') ==7)
    mkdir(dem_save_root);
end
dem_save = sprintf('%sdem_%04d_yearly.mat',dem_save_root,iyear)';

% read cropland fraction
crop_fraction_fn = sprintf('%sdata/Global_cropland_3km_2000_2019.mat',root_dir);
crop_fraction = loadMatData(crop_fraction_fn);

crop_mask = lc==12 | lc==14 | crop_fraction>=40;
non_crop_mask = lc>0 & lc ~= 13 & ~crop_mask;  % exclude water and crop and urban

path_in = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_input_eramodis/yearly/diurnal/');
path_out = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_input_eramodis/yearly/diurnal_nearbyall/');
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

%% SH & LE
[sshf_save,sshf]=load_data(path_in,path_out,iyear,'sshf');
[slhf_save,slhf]=load_data(path_in,path_out,iyear,'slhf');

%% emis
[EMIS_save,EMIS]=load_data(path_in,path_out,iyear,'EMIS');

%% 2t skt
[t2_save,t2]=load_data(path_in,path_out,iyear,'2t');
[skt_save,skt]=load_data(path_in,path_out,iyear,'skt');

%%  6 & 7  & 8 Shortwave in,  Longwave in, Albedo,
[ssrd_save,ssrd]=load_data(path_in,path_out,iyear,'ssrd');
[strd_save,strd]=load_data(path_in,path_out,iyear,'strd');
[fal_save,fal]=load_data(path_in,path_out,iyear,'fal');

%% 9 specific humidity
[Qair_save,Qair]=load_data(path_in,path_out,iyear,'Qair');

%% 10 & 11  air density, air pressure
[RHOA_save,RHOA]=load_data(path_in,path_out,iyear,'RHOA');
[sp_save,sp]=load_data(path_in,path_out,iyear,'sp');

%% 12 ground heat flux
[GH_save,GH]=load_data(path_in,path_out,iyear,'GH');

%% 13 & 14 aerodynamic resistance, surface resistance
[RA_save,RA]=load_data(path_in,path_out,iyear,'RA');
[RS_save,RS]=load_data(path_in,path_out,iyear,'RS');



%% common mask
common_mask = ~isnan(sshf) & ~isnan(slhf) & ~isnan(EMIS) & ~isnan(t2) & ~isnan(skt) & ~isnan(ssrd) & ~isnan(strd)...
    & ~isnan(fal) & ~isnan(Qair) & ~isnan(RHOA) & ~isnan(sp) & ~isnan(GH) & ~isnan(RA) & ~isnan(RS);

sshf(~common_mask) = nan;
slhf(~common_mask) = nan;
EMIS(~common_mask) = nan;
t2(~common_mask) = nan;
skt(~common_mask) = nan;
ssrd(~common_mask) = nan;
strd(~common_mask) = nan;
fal(~common_mask) = nan;
Qair(~common_mask) = nan;
RHOA(~common_mask) = nan;
sp(~common_mask) = nan;
GH(~common_mask) = nan;
RA(~common_mask) = nan;
RS(~common_mask) = nan;

dem(~common_mask)=nan;

%% initialize
sshf_new = nan(ipt.row,ipt.col,'single');
slhf_new = nan(ipt.row,ipt.col,'single');
EMIS_new = nan(ipt.row,ipt.col,'single');
t2_new =  nan(ipt.row,ipt.col,'single');
skt_new =  nan(ipt.row,ipt.col,'single');
ssrd_new =  nan(ipt.row,ipt.col,'single');
strd_new =  nan(ipt.row,ipt.col,'single');
fal_new =  nan(ipt.row,ipt.col,'single');
Qair_new =  nan(ipt.row,ipt.col,'single');
RHOA_new =  nan(ipt.row,ipt.col,'single');
sp_new =  nan(ipt.row,ipt.col,'single');
GH_new =  nan(ipt.row,ipt.col,'single');
RA_new =  nan(ipt.row,ipt.col,'single');
RS_new =  nan(ipt.row,ipt.col,'single');
dem_new = nan(ipt.row,ipt.col,'single');


pixel2search = find(srcwin>0 & common_mask);
n_total_pixel = numel(pixel2search);
icount=0;
for id_pixel = pixel2search'
    icount = icount+1;
    [irow,icol]=ind2sub([ipt.row,ipt.col],id_pixel);
       
    offset = (srcwin(irow,icol)-1)/2;
    if irow-offset<1
        irow_start=1;
    else
        irow_start = irow-offset;
    end
    
    if irow+offset>ipt.row
        irow_end = ipt.row;
    else
        irow_end = irow+offset;
    end
    
    if icol-offset<1
        icol_start=1;
    else
        icol_start = icol-offset;
    end
    
    if icol+offset>ipt.col
        icol_end = ipt.col;
    else
        icol_end = icol+offset;
    end
       
    roi_mask = non_crop_mask(irow_start:irow_end,icol_start:icol_end); % data with value 1 in this mask will be averaged
    
    crop_dem = dem(irow,icol);
    dem_roi = dem(irow-offset:irow+offset,icol-offset:icol+offset);
    dem_roi(~roi_mask)=nan;
    delta_dem = crop_dem - dem_roi;
    roi_mask(abs(delta_dem)>100)=false;           
        
    % SH LE
    sshf_roi = sshf(irow_start:irow_end,icol_start:icol_end);
    sshf_new(irow,icol)=mean(sshf_roi(roi_mask),'omitnan');
    
    slhf_roi = slhf(irow_start:irow_end,icol_start:icol_end);
    slhf_new(irow,icol)=mean(slhf_roi(roi_mask),'omitnan');
    
    %EMIS
    EMIS_roi = EMIS(irow_start:irow_end,icol_start:icol_end);
    EMIS_new(irow,icol)=mean(EMIS_roi(roi_mask),'omitnan');
    
    % 2t
    t2_roi = t2(irow_start:irow_end,icol_start:icol_end);
    t2_pixel = t2_roi(roi_mask);
    npixel = sum(~isnan(t2_pixel));
    t2_new(irow,icol)=(sum(t2_pixel.^4,'omitnan')/npixel)^0.25;
    
    % skt
    skt_roi = skt(irow_start:irow_end,icol_start:icol_end);
    skt_pixel = skt_roi(roi_mask);
    EMIS_pixel = EMIS_roi(roi_mask);
    
    temp=EMIS_pixel.*skt_pixel.^4;
    npixel = sum(~isnan(temp));
    skt_new(irow,icol)=(sum(temp,'omitnan')/npixel/ EMIS_new(irow,icol))^0.25;
    
    %ssrd
    ssrd_roi = ssrd(irow_start:irow_end,icol_start:icol_end);
    ssrd_new(irow,icol)=mean(ssrd_roi(roi_mask),'omitnan');
    
    %strd
    strd_roi = strd(irow_start:irow_end,icol_start:icol_end);
    strd_new(irow,icol)=mean(strd_roi(roi_mask),'omitnan');
    
    %fal
    fal_roi = fal(irow_start:irow_end,icol_start:icol_end);
    fal_new(irow,icol)=mean(fal_roi(roi_mask),'omitnan');
    
    %Qair
    Qair_roi = Qair(irow_start:irow_end,icol_start:icol_end);
    Qair_new(irow,icol)=mean(Qair_roi(roi_mask),'omitnan');
    
    %sp
    sp_roi = sp(irow_start:irow_end,icol_start:icol_end);
    sp_new(irow,icol)=mean(sp_roi(roi_mask),'omitnan');
    
    %RHOA
    RHOA_roi = RHOA(irow_start:irow_end,icol_start:icol_end);
    RHOA_new(irow,icol)=mean(RHOA_roi(roi_mask),'omitnan');
    
    %GH
    GH_roi = GH(irow_start:irow_end,icol_start:icol_end);
    GH_new(irow,icol)=mean(GH_roi(roi_mask),'omitnan');
    
    %RA
    RA_roi = RA(irow_start:irow_end,icol_start:icol_end);
    RA_new(irow,icol)=mean(RA_roi(roi_mask),'omitnan');
    
    %RS
    RS_roi = RS(irow_start:irow_end,icol_start:icol_end);
    RS_new(irow,icol)=mean(RS_roi(roi_mask),'omitnan');
    
    % dem
    dem_roi = dem(irow_start:irow_end,icol_start:icol_end);
    dem_new(irow,icol)=mean(dem_roi(roi_mask),'omitnan');
    
    pct_finished = icount / n_total_pixel * 100;
    if mod(pct_finished,5)==0
        fprintf('finished %f%%\n',pct_finished)
    end
end


%% SH & LE
save(sshf_save,'sshf_new','-v7.3');
save(slhf_save,'slhf_new','-v7.3');

%% emis
save(EMIS_save,'EMIS_new','-v7.3');

%% 2t skt
save(t2_save,'t2_new','-v7.3');
save(skt_save,'skt_new','-v7.3');

%%  6 & 7  & 8 Shortwave in,  Longwave in, Albedo,
save(ssrd_save,'ssrd_new','-v7.3');
save(strd_save,'strd_new','-v7.3');
save(fal_save,'fal_new','-v7.3');

%% 9 specific humidity
save(Qair_save,'Qair_new','-v7.3');

%% 10 & 11  air density, air pressure
save(RHOA_save,'RHOA_new','-v7.3');
save(sp_save,'sp_new','-v7.3');

%% 12 ground heat flux
save(GH_save,'GH_new','-v7.3');

%% 13 & 14 aerodynamic resistance, surface resistance
save(RA_save,'RA_new','-v7.3');
save(RS_save,'RS_new','-v7.3');

%% 15 dem
save(dem_save,'dem_new','-v7.3');

toc
end


%% load data
function [save_path,data]=load_data(path_in,path_out,iyear,varname)
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
save_root = sprintf('%s%s/',path_out,varname);
if ~(exist(save_root,'dir') ==7)
    mkdir(save_root);
end
save_path = sprintf('%s%s/%s_%04d_yearly.mat',path_out,varname,varname,iyear);
data = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',path_in,varname,varname,iyear));
end
