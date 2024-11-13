function A000_main_function()

%% data processing
A001_lc_analysis_c61()

A002_get_cropland_adjacent()

A003_get_dominant_surrounding_biome()

for iyear = 2001:2023
    A004_prepare_TRM_inputs_yearly_diurnal_MODIS_ERA(iyear)
    A005_prepare_TRM_inputs_yearly_diurnal_MODIS_ERA_adj(iyear)
    A006_Get_TRM_outputs(iyear)   
end

end

