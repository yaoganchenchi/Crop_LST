function A001_lc_analysis_c61() % c6 lc
%% this script will extract the mode land cover based on 21 MODIS land cover maps

ipt.tool_path = '/projectsp/f_cc2127_1/chenchi/matlab_tool/';
addpath(ipt.tool_path);
ipt.hdf_path = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/data/'; % path to your MCD12C1 data
ipt.mat_path = '/projectsp/f_cc2127_1/chenchi/Data_Archive/MCD12C1/output_c61/';
ipt.Nrow = 3600;
ipt.Ncol = 7200;
ipt.syear=2001;
ipt.eyear=2005;
ipt.nlayer = ipt.eyear-ipt.syear+1;

get_mode_lc(ipt)

end

%% get MODE LC
function get_mode_lc(ipt)
addpath(ipt.tool_path);
data = zeros(ipt.Nrow, ipt.Ncol,ipt.nlayer,'uint8');
ict = 0;
for iyear = ipt.syear:ipt.eyear
    ict= ict+1;
    fn = sprintf('%sMCD12C1.A%04d001.061.*',ipt.hdf_path,iyear);
    file =dir (fn);
    fn = sprintf('%s%s',ipt.hdf_path,file(1).name);
    temp = hdfread(fn, '/MOD12C1/Data Fields/Majority_Land_Cover_Type_1', 'Index', {[1  1],[1  1],[3600  7200]});
    data(:,:,ict) = temp;

end
data=single(data);
data(data>16)=nan;
mode_lc = mode(data,3);

savename = sprintf('%sMCD12C1.CMG005.C61.IGBP.MODE.LC.%04d.%04d.mat',ipt.mat_path,ipt.syear,ipt.eyear);
save(savename,'mode_lc','-v7.3')

end