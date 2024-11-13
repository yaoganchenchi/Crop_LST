function A002_get_cropland_adjacent()
%% this script will determine the search window size based on MODIS land cover
root_dir = '/projectsp/f_cc2127_1/chenchi/LCLUC-crop/';
addpath(sprintf('%scode/matlab_tools/',root_dir));

%% read modis lc
fn_lc = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/output_c61/MCD12C1.CMG005.C61.IGBP.MODE.LC.2001.2021.mat';
lc = loadMatData(fn_lc);

%% read cropland fraction
output_dir  = sprintf('%sdata/',root_dir);
crop_fraction_fn = sprintf('%sGlobal_cropland_3km_2000_2019.mat',output_dir);
crop_fraction = loadMatData(crop_fraction_fn);

dem_fn = sprintf('%sdata/GMTED2010_DEM_CMG005.tif',root_dir);
dem = readgeoraster(dem_fn);

%% 
crop_mask = lc==12 | lc==14 | crop_fraction>=40;
non_crop_mask = lc>0 & lc ~= 13 & ~crop_mask;  % exclude water and crop and urban

%% determine moving window size by MCD12
srch_win=nan(3600,7200);
for irow=1:3600
    for icol=1:7200
        if ~crop_mask(irow,icol)
            continue
        end

        %% does not check boundary, maybe unnecessary
        min_size =3; % base
        while isnan(srch_win(irow,icol))
            min_size=min_size+2;  % first window will be 3
            offset = (min_size -1)/2;
            roi = non_crop_mask(irow-offset:irow+offset,icol-offset:icol+offset);
            crop_dem = dem(irow,icol);
            
            dem_roi = dem(irow-offset:irow+offset,icol-offset:icol+offset);
            dem_roi(~roi)=nan;
            delta_dem = crop_dem - dem_roi;
            roi(abs(delta_dem)>100)=false;           
            n_non_veg = sum(roi(:));
            if  n_non_veg >0
                srch_win(irow,icol)=min_size;    
            end
        end
    end
end
savename = sprintf('%ssearch_window_size_MCD12.mat',output_dir);
save(savename,'srch_win','-v7.3')

end



