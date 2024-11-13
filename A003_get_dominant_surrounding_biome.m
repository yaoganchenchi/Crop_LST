function A003_get_dominant_surrounding_biome()
%% this script returns the dominant surrounding biome of each cropland pixel (with a mix of multilpe biomes)

root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));
addpath(sprintf('%scode/TRM/',root_dir));

fn_lc = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/output_c61/MCD12C1.CMG005.C61.IGBP.MODE.LC.2001.2021.mat';
fn_srcwin = sprintf('%ssearch_window_size_MCD12.mat','/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/');

lc = loadMatData(fn_lc);
srcwin = loadMatData(fn_srcwin);

fn_dir  = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/');
crop_fraction_fn = sprintf('%sGlobal_cropland_3km_2000_2019.mat',fn_dir);
crop_fraction = loadMatData(crop_fraction_fn);

crop_mask = lc==12 | lc==14 | crop_fraction>=40;
non_crop_mask = lc>0 & lc ~= 13 & ~crop_mask;  % exclude water and crop

dem_fn = sprintf('%sdata/GMTED2010_DEM_CMG005.tif',root_dir);
dem = readgeoraster(dem_fn);

lc(~non_crop_mask)=nan;
lc_broad = nan(size(lc));

% define the broad biome types
lc_broad((lc>=1 & lc<=5) | (lc>=8 & lc<=9))=1; % trees
lc_broad(lc==10 | lc==11)=2; % Grass & Wetlands
lc_broad((lc>=6 & lc<=7)| (lc==15 | lc ==16))=3; % Shrubs, barren, ice, no urban

pixel2search = find(srcwin>0);

path_out = sprintf('/projectsp/f_cc2127_1/chenchi/LCLUC-crop/data/');
if ~(exist(path_out,'dir') ==7)
    mkdir(path_out);
end

ipt.R = georefcells([-90 90],[-180 180],[1800,3600],'ColumnsStartFrom','north');
ipt.row=3600;
ipt.col=7200;

surrounding_lc =  nan(ipt.row,ipt.col,'single');

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
    
    roi = lc_broad(irow_start:irow_end,icol_start:icol_end); % data with value 1 in this mask will be averaged
    
    crop_dem = dem(irow,icol);    
    dem_roi = dem(irow-offset:irow+offset,icol-offset:icol+offset);
    
    dem_roi(isnan(roi))=nan;
    delta_dem = crop_dem - dem_roi;
    roi(abs(delta_dem)>100)=nan;
       
    surrounding_lc(irow,icol) = mode(roi(:));

    pct_finished = icount / n_total_pixel * 100;
    if mod(pct_finished,5)==0
        fprintf('finished %f%%\n',pct_finished)
    end
end

save_path = sprintf('%sdominant_surrounding_lc.mat',path_out);
save(save_path,'surrounding_lc','-v7.3');

end