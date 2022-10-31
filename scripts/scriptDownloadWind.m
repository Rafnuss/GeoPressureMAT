
%% Load data

project="StudyStarling";
load("../data/processedData"+project+".mat")


addpath(genpath('../functions'))
pyenv('Version','/usr/local/opt/python@3.8/bin/python3.8')



%% Download

c = py.cdsapi.Client();
possible_pressure=[1 2 3 5 7 10 20 30 50 70 100:25:250 300:50:750 775:25:1000];

for lt=1:height(tblLog)
    mkdir("../data/ECMWF/" + tblLog.GDL_ID{lt})
    for i_s=1:height(sta{lt})-1
        
        t = dateshift(sta{lt}.end(i_s),'start','hour')-1/24:1/24:dateshift(sta{lt}.start(i_s+1),'end','hour')+1/24;

        id_pres = t(1)<raw{lt}.pressure.date & raw{lt}.pressure.date < t(end);
        [~,id_min] = find(min(raw{lt}.pressure.obsWithOutliars(id_pres))>possible_pressure,1,'last');
        id_min = min(id_min,numel(possible_pressure)-1);
        if any(max(raw{lt}.pressure.obs(id_pres))<possible_pressure)
            [~,id_max] = find(max(raw{lt}.pressure.obs(id_pres))<possible_pressure,1);
        else
            id_max=numel(possible_pressure);
        end
        id_prespos=id_min:id_max;

        if numel(id_prespos)==0
            keyboard
        end
        
        a="/usr/local/opt/python@3.8/bin/python3.8 -c """...
            +"import cdsapi;"...
            +"c = cdsapi.Client();"...
            +"data = c.retrieve('reanalysis-era5-pressure-levels',"...
            +"{"...
            +"    'product_type': 'reanalysis',"...
            +"    'format': 'netcdf',"...
            +"    'variable': ['u_component_of_wind', 'v_component_of_wind'],"...
            +"    'year': ["+replace(num2str(sort(unique(year(t)))),'  ',',')+"],"...
            +"    'month': ["+replace(num2str(sort(unique(month(t)))),'  ',',')+"],"...
            +"    'day': ["+replace(num2str(sort(unique(day(t)))),'  ',',')+"],"...
            +"    'time': ["+replace(num2str(sort(unique(hour(t)))),'  ',',')+"],"...
            +"    'pressure_level': ["+replace(num2str(possible_pressure(id_prespos)),'  ',',')+"],"...
            +"    'area': ["+tblLog.bndy_N(lt)+", "+tblLog.bndy_W(lt)+", "+tblLog.bndy_S(lt)+", "+tblLog.bndy_E(lt)+"],"...
            +"},"...
            +"'../data/ECMWF/" + tblLog.GDL_ID{lt} + "/" + "sta_" + num2str(i_s) +".nc')"...
            +""" &";
        system(a)
    end
end





%% Group migration period for fewer files
% explore optimal threashold of migration group duration
% tmp = cellfun(@(x) datenum(diff(x.end(1:end-1)))',sta,'UniformOutput',false);
% figure; hold on; histogram([tmp{:}],1:15)
% thr=1;


%% Download
% Version grouping migratory period

% c = py.cdsapi.Client();
% possible_pressure=[1 2 3 5 7 10 20 30 50 70 100:25:250 300:50:750 775:25:1000];
% 
% for lt=1:height(tblLog)
%     mkdir("ECMWF/" + tblLog.GDL_ID{lt})
% 
%     % Compute the start and end of each migration bout
%     migStart = dateshift(sta{lt}.end(1:end-1),'start','hour')-1/24;
%     migEnd = dateshift(sta{lt}.start(2:end),'end','hour')+1/24;
% 
%     % Group migration seperated by thr days at least
%     T = clusterdata(datenum(migStart),'Criterion','distance','Cutoff',thr);
%     % figure; plot(migStart,T,'o')
% 
%     parfor i_t=2:numel(unique(T))
% 
%         t = min(migStart(i_t==T)):1/24:max(migEnd(i_t==T));
% 
%         filename =  datestr(t(1),'yyyymmddHH') + "-" + datestr(t(end),'yyyymmddHH') +".nc";
% 
%         id_pres = t(1)<raw{lt}.pressure.date & raw{lt}.pressure.date < t(end);
%         id_prespos=min(raw{lt}.pressure.obs(id_pres))<possible_pressure & max(raw{lt}.pressure.obs(id_pres))>possible_pressure;
%         id_prespos=movmean(id_prespos,3)>0; % select the pressure before and after in addition.
% 
%         
%         system("/usr/local/opt/python@3.8/bin/python3.8 -c """...
%             +"import cdsapi;"...
%             +"c = cdsapi.Client();"...
%             +"data = c.retrieve('reanalysis-era5-pressure-levels',"...
%             +"{"...
%             +"    'product_type': 'reanalysis',"...
%             +"    'format': 'netcdf',"...
%             +"    'variable': ['u_component_of_wind', 'v_component_of_wind'],"...
%             +"    'year': ["+replace(num2str(sort(unique(year(t)))),'  ',',')+"],"...
%             +"    'month': ["+replace(num2str(sort(unique(month(t)))),'  ',',')+"],"...
%             +"    'day': ["+replace(num2str(sort(unique(day(t)))),'  ',',')+"],"...
%             +"    'time': ["+replace(num2str(sort(unique(hour(t)))),'  ',',')+"],"...
%             +"    'pressure_level': ["+replace(num2str(possible_pressure(id_prespos)),'  ',',')+"],"...
%             +"    'area': ["+tblLog.bndy_N(lt)+", "+tblLog.bndy_W(lt)+", "+tblLog.bndy_S(lt)+", "+tblLog.bndy_E(lt)+"],"...
%             +"},"...
%             +"'ECMWF/" + tblLog.GDL_ID{lt} + "/" + filename + "')"...
%             +"""")
%     end
% end