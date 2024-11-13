function A006_Get_TRM_outputs(iyear) % diurnal
if ~exist('iyear','var')
    iyear = 2022;
end

tic

time_scale1 = 'yearly';
fprintf('iyear=%04d, season: %s\n',iyear,time_scale1);

% by year, and by month
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
addpath(sprintf('%scode/TRM/',root_dir));

%% yearly, diurnal
path_in = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_input_eramodis/');

% surrounding veg, get dTs to dbiophysical
time_scale2 = 'diurnal_nearbyall';  %time scale 2
dir0 = sprintf('%s%s/%s/',path_in,time_scale1,time_scale2);
ipt1.Ts = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'skt','skt',iyear));
ipt1.SWin = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'ssrd','ssrd',iyear));
ipt1.LWin = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'strd','strd',iyear));
ipt1.Ta   = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'2t','2t',iyear));
ipt1.qa   = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'Qair','Qair',iyear));
ipt1.albedo = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'fal','fal',iyear));
ipt1.ra = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RA','RA',iyear));
ipt1.rs = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RS','RS',iyear));
ipt1.emis = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'EMIS','EMIS',iyear));
ipt1.rhoa = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RHOA','RHOA',iyear));
ipt1.Ps = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'sp','sp',iyear));
ipt1.G  = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'GH','GH',iyear));
ipt1.SH = -loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'sshf','sshf',iyear));
ipt1.LE = -loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'slhf','slhf',iyear));
ipt1.dem = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'dem','dem',iyear));

mask1 = ~isnan(ipt1.SWin) & ~isnan(ipt1.LWin) & ~isnan(ipt1.Ta) & ~isnan(ipt1.qa)...
& ~isnan(ipt1.albedo) & ~isnan(ipt1.ra) & ~isnan(ipt1.rs) & ~isnan(ipt1.emis)...
& ~isnan(ipt1.rhoa) & ~isnan(ipt1.Ps) & ~isnan(ipt1.G) & ~isnan(ipt1.Ts);

% crop , get dTs to dbiophysical
time_scale2 = 'diurnal';  %time scale 2
dir0 = sprintf('%s%s/%s/',path_in,time_scale1,time_scale2);
path_out = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/TRM_output/%s/%s/%04d/',time_scale1,time_scale2,iyear);
if ~(exist(path_out,'dir') ==7)
    mkdir(path_out);
end
ipt2.Ts = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'skt','skt',iyear));

fn_lc = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/output_c61/MCD12C1.CMG005.C61.IGBP.MODE.LC.2001.2021.mat';
lc = loadMatData(fn_lc);
fn_dir  = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/');
crop_fraction_fn = sprintf('%sGlobal_cropland_3km_2000_2019.mat',fn_dir);
crop_fraction = loadMatData(crop_fraction_fn);

crop_mask = lc==12 | lc==14 | crop_fraction>=40;  %  crop
ipt2.Ts(~crop_mask)=nan;

ipt2.SWin = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'ssrd','ssrd',iyear));
ipt2.LWin = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'strd','strd',iyear));
ipt2.Ta  = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'2t','2t',iyear));
ipt2.qa    = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'Qair','Qair',iyear));
ipt2.albedo = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'fal','fal',iyear));
ipt2.ra = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RA','RA',iyear));
ipt2.rs = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RS','RS',iyear));
ipt2.emis = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'EMIS','EMIS',iyear));
ipt2.rhoa = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'RHOA','RHOA',iyear));
ipt2.Ps = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'sp','sp',iyear));
ipt2.G  = loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'GH','GH',iyear));
ipt2.SH = -loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'sshf','sshf',iyear));
ipt2.LE = -loadMatData(sprintf('%s%s/%s_%04d_yearly.mat',dir0,'slhf','slhf',iyear));

dem_fn = sprintf('%sdata/GMTED2010_DEM_CMG005.tif',root_dir);
ipt2.dem = readgeoraster(dem_fn);

mask2 = ~isnan(ipt2.SWin) & ~isnan(ipt2.LWin) & ~isnan(ipt2.Ta) & ~isnan(ipt2.qa)...
& ~isnan(ipt2.albedo) & ~isnan(ipt2.ra) & ~isnan(ipt2.rs) & ~isnan(ipt2.emis)...
& ~isnan(ipt2.rhoa) & ~isnan(ipt2.Ps) & ~isnan(ipt2.G) & ~isnan(ipt2.Ts);

% data filtering
final_mask = mask1 & mask2 & crop_mask;

ipt1.Ts(~final_mask)=nan;
ipt1.SWin(~final_mask)=nan;
ipt1.LWin(~final_mask)=nan;
ipt1.Ta(~final_mask)=nan;
ipt1.qa(~final_mask)=nan;
ipt1.albedo(~final_mask)=nan;
ipt1.ra(~final_mask)=nan;
ipt1.rs(~final_mask)=nan;
ipt1.emis(~final_mask)=nan;
ipt1.rhoa(~final_mask)=nan;
ipt1.Ps(~final_mask)=nan;
ipt1.G(~final_mask)=nan;
ipt1.dem(~final_mask)=nan;

ipt2.Ts(~final_mask)=nan;
ipt2.SWin(~final_mask)=nan;
ipt2.LWin(~final_mask)=nan;
ipt2.Ta(~final_mask)=nan;
ipt2.qa(~final_mask)=nan;
ipt2.albedo(~final_mask)=nan;
ipt2.ra(~final_mask)=nan;
ipt2.rs(~final_mask)=nan;
ipt2.emis(~final_mask)=nan;
ipt2.rhoa(~final_mask)=nan;
ipt2.Ps(~final_mask)=nan;
ipt2.G(~final_mask)=nan;
ipt2.dem(~final_mask)=nan;

[sen1] = TRM(ipt1.SWin, ipt1.LWin, ipt1.Ta, ipt1.qa, ipt1.albedo, ipt1.ra, ipt1.rs, ipt1.emis, ipt1.rhoa, ipt1.Ps, ipt1.G);
[sen2] = TRM(ipt2.SWin, ipt2.LWin, ipt2.Ta, ipt2.qa, ipt2.albedo, ipt2.ra, ipt2.rs, ipt2.emis, ipt2.rhoa, ipt2.Ps, ipt2.G);

%% optimize m
block_size = 200;
[nrow,ncol] = size(ipt1.Ts);
true_delta_Ts = (ipt2.Ts-ipt1.Ts);

delta_var.albedo = (ipt2.albedo-ipt1.albedo);
delta_var.ra = (ipt2.ra-ipt1.ra);
delta_var.rs = (ipt2.rs-ipt1.rs);
delta_var.emis = (ipt2.emis-ipt1.emis);
delta_var.SWin = (ipt2.SWin-ipt1.SWin);
delta_var.LWin = (ipt2.LWin-ipt1.LWin);
delta_var.qa = (ipt2.qa-ipt1.qa);
delta_var.Ta = (ipt2.Ta-ipt1.Ta);
delta_var.G = (ipt2.G-ipt1.G);

delta_var.SH = (ipt2.SH-ipt1.SH);
delta_var.LE = (ipt2.LE-ipt1.LE);
delta_var.LAI=(ipt2.LAI-ipt1.LAI);
delta_var.uv10 = (ipt2.u10.^2 + ipt2.v10.^2).^0.5 - (ipt1.u10.^2 + ipt1.v10.^2).^0.5;
delta_var.eraskt = (ipt2.eraskt-ipt1.eraskt);
delta_var.dem = (ipt2.dem-ipt1.dem);

fields_delta = fieldnames(delta_var);
fields_sen = fieldnames(sen1);

n_row_block = nrow/block_size;
n_col_block = ncol/block_size;

m_opt = nan(nrow,ncol);
for i_col_block = 1:n_col_block
    fprintf('processing iblock %d/36\n',i_col_block)
    for i_row_block = 1:n_row_block
        row_idx = (i_row_block-1)*block_size+1 : (i_row_block)*block_size;
        col_idx = (i_col_block-1)*block_size+1 : (i_col_block)*block_size;
        valid_roi = final_mask(row_idx,col_idx);
        
        % get ROI for true delta Ts
        true_delta_Ts_roi = true_delta_Ts(row_idx,col_idx);
        % get ROI for delta
        for ii = 1: length(fields_delta)
            this_field = fields_delta{ii};
            delta_var_roi.(this_field) = delta_var.(this_field)(row_idx,col_idx);
        end
        % get ROI for sen
        for ii = 1: length(fields_sen)
            this_field = fields_sen{ii};
            sen1_roi.(this_field) = sen1.(this_field)(row_idx,col_idx);
            sen2_roi.(this_field) = sen2.(this_field)(row_idx,col_idx);
        end

        m = 0:0.1:10;
        Residual = nan(block_size,block_size,length(m));
        for i = 1:length(m)
            sen_roi.dTs_dALBEDO = (sen1_roi.dTs_dALBEDO + m(i).*sen2_roi.dTs_dALBEDO)./(1+m(i));
            sen_roi.dTs_dRA = (sen1_roi.dTs_dRA + m(i).*sen2_roi.dTs_dRA)./(1+m(i));
            sen_roi.dTs_dRS = (sen1_roi.dTs_dRS + m(i).*sen2_roi.dTs_dRS)./(1+m(i));
            sen_roi.dTs_dEMIS = (sen1_roi.dTs_dEMIS + m(i).*sen2_roi.dTs_dEMIS)./(1+m(i));
            sen_roi.dTs_dSWin = (sen1_roi.dTs_dSWin + m(i).*sen2_roi.dTs_dSWin)./(1+m(i));
            sen_roi.dTs_dLWin = (sen1_roi.dTs_dLWin + m(i).*sen2_roi.dTs_dLWin)./(1+m(i));
            sen_roi.dTs_dQA = (sen1_roi.dTs_dQA + m(i).*sen2_roi.dTs_dQA)./(1+m(i));
            sen_roi.dTs_dTA = (sen1_roi.dTs_dTA + m(i).*sen2_roi.dTs_dTA)./(1+m(i));
            sen_roi.dTs_dG = (sen1_roi.dTs_dG + m(i).*sen2_roi.dTs_dG)./(1+m(i));

            delta_Ts_calculated_roi = sen_roi.dTs_dALBEDO.*delta_var_roi.albedo + sen_roi.dTs_dRA.*delta_var_roi.ra + sen_roi.dTs_dRS.*delta_var_roi.rs +...
                sen_roi.dTs_dEMIS.*delta_var_roi.emis + sen_roi.dTs_dSWin.*delta_var_roi.SWin + sen_roi.dTs_dLWin.*delta_var_roi.LWin +...
                sen_roi.dTs_dQA.*delta_var_roi.qa + sen_roi.dTs_dTA.* delta_var_roi.Ta + sen_roi.dTs_dG.*delta_var_roi.G;
            
            Residual(:,:,i) = abs(delta_Ts_calculated_roi-true_delta_Ts_roi);
        end
        [Residual_min, Residual_index] = min(Residual,[],3);
        m_opt_roi = m(Residual_index);
        m_opt_roi(~valid_roi)=nan;
        m_opt(row_idx,col_idx) = m_opt_roi;
    end
end

sen_opt.dTs_dALBEDO = (sen1.dTs_dALBEDO + m_opt.*sen2.dTs_dALBEDO)./(1+m_opt);
sen_opt.dTs_dRA = (sen1.dTs_dRA + m_opt.*sen2.dTs_dRA)./(1+m_opt);
sen_opt.dTs_dRS = (sen1.dTs_dRS + m_opt.*sen2.dTs_dRS)./(1+m_opt);
sen_opt.dTs_dEMIS = (sen1.dTs_dEMIS + m_opt.*sen2.dTs_dEMIS)./(1+m_opt);
sen_opt.dTs_dSWin = (sen1.dTs_dSWin + m_opt.*sen2.dTs_dSWin)./(1+m_opt);
sen_opt.dTs_dLWin = (sen1.dTs_dLWin + m_opt.*sen2.dTs_dLWin)./(1+m_opt);
sen_opt.dTs_dQA = (sen1.dTs_dQA + m_opt.*sen2.dTs_dQA)./(1+m_opt);
sen_opt.dTs_dTA = (sen1.dTs_dTA + m_opt.*sen2.dTs_dTA)./(1+m_opt);
sen_opt.dTs_dG = (sen1.dTs_dG + m_opt.*sen2.dTs_dG)./(1+m_opt);

% save 1
output.TRM_m_sp1_delta_Ts = sen_opt.dTs_dALBEDO.*delta_var.albedo + sen_opt.dTs_dRA.*delta_var.ra + sen_opt.dTs_dRS.*delta_var.rs +...
    sen_opt.dTs_dEMIS.*delta_var.emis + sen_opt.dTs_dSWin.*delta_var.SWin + sen_opt.dTs_dLWin.*delta_var.LWin +...
    sen_opt.dTs_dQA.*delta_var.qa + sen_opt.dTs_dTA.* delta_var.Ta + sen_opt.dTs_dG.*delta_var.G;
output.Ts_crop_MODIS_true = ipt2.Ts;
output.Ts_other_MODIS_true = ipt1.Ts;
output.true_delta_Ts = true_delta_Ts;

output.bias_TRM_m_sp1 = (output.TRM_m_sp1_delta_Ts - true_delta_Ts);  % k  %/true_delta_Ts;
output.bias_TRM_m_sp1_rel = (output.TRM_m_sp1_delta_Ts - true_delta_Ts) ./true_delta_Ts * 100;  % ;

output.Ts_crop_TRM = sen2.Ts;
output.Ts_other_TRM = sen1.Ts;
output.TRM_direct_delta_Ts = sen2.Ts-sen1.Ts;

savename = sprintf('%sdelta_Ts_%04d.mat',path_out,iyear);
save(savename,'output','-v7.3');

% save 2
breakdown.dTs_albedo = sen_opt.dTs_dALBEDO.*delta_var.albedo;
breakdown.dTs_RA = sen_opt.dTs_dRA.*delta_var.ra;
breakdown.dTs_RS = sen_opt.dTs_dRS.*delta_var.rs;
breakdown.dTs_EMIS = sen_opt.dTs_dEMIS.*delta_var.emis;
breakdown.dTs_SWin = sen_opt.dTs_dSWin.*delta_var.SWin;
breakdown.dTs_LWin = sen_opt.dTs_dLWin.*delta_var.LWin;
breakdown.dTs_QA = sen_opt.dTs_dQA.*delta_var.qa;
breakdown.dTs_TA = sen_opt.dTs_dTA.*delta_var.Ta;
breakdown.dTs_G = sen_opt.dTs_dG.*delta_var.G;
savename = sprintf('%sdelta_Ts_breakdown_%04d.mat',path_out,iyear);
save(savename,'breakdown','-v7.3');

% save 3
savename = sprintf('%sdelta_var_%04d.mat',path_out,iyear);
save(savename,'delta_var','-v7.3');

% save 4
out1.m_opt =m_opt;
out1.sen_opt=sen_opt;
savename = sprintf('%ssen_opt_m_%04d.mat',path_out,iyear);
save(savename,'out1','-v7.3');

%save 5
out2.sen_crop = sen2;
out2.sen_other = sen1;
savename = sprintf('%ssen_original_%04d.mat',path_out,iyear);
save(savename,'out2','-v7.3');

toc

end

