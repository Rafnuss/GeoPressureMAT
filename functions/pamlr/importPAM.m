function s = importPAM(filepath,varargin)
%importPAM Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(varargin{1}) || isnat(varargin{1})
    s_d = datetime('2000-01-01');
else
    s_d = varargin{1};
end

if nargin<2 || isempty(varargin{2}) || isnat(varargin{2})
    e_d = datetime('2100-01-01');
else
    e_d = varargin{2};
end

tmp = strsplit(filepath,'/');
filename = [filepath '/' tmp{end}];

% Read settings
s = jsondecode(fileread([filename '.settings']));
s.StarttimeRTC = datetime(s.StarttimeRTC,'InputFormat','dd.MM.yyyy HH:mm:ss');
s.LoggingStart = datetime(s.LoggingStart,'InputFormat','dd.MM.yyyy');
s.StopTimeRTC = datetime(s.StopTimeRTC,'InputFormat','dd.MM.yyyy HH:mm:ss');
s.StopTimeReference = datetime(s.StopTimeReference,'InputFormat','dd.MM.yyyy HH:mm:ss');
s.StopTimeRTC = datetime(s.StopTimeRTC,'InputFormat','dd.MM.yyyy HH:mm:ss');

% Reading light
dataArray = readFile([filename '.glf'],'%17{dd.MM.yyyy HH:mm}D %f %f %f %f',s_d,e_d);
s.light.date = dataArray{1};
s.light.obs = dataArray{2};

% Reading Pressure
dataArray = readFile([filename '.pressure'],'%17{dd.MM.yyyy HH:mm}D %f',s_d,e_d);
s.pressure.date = dataArray{1};
s.pressure.obs = dataArray{2};

% Reading Acceleration
dataArray = readFile([filename '.acceleration'],'%17{dd.MM.yyyy HH:mm}D %f %f',s_d,e_d);
s.acceleration.date = dataArray{1};
s.acceleration.pit = dataArray{2};
s.acceleration.act = dataArray{3};

% Reading temperature
if exist([filename '.AirTemperature'],'file')
    dataArray = readFile([filename '.AirTemperature'],'%17{dd.MM.yyyy HH:mm}D %f',s_d,e_d);
else
    dataArray = readFile([filename '.temperature'],'%17{dd.MM.yyyy HH:mm}D %f',s_d,e_d);
end
s.temperature.date = dataArray{1};
s.temperature.obs = dataArray{2};

% Reading magnetic
dataArray = readFile([filename '.magnetic'],'%17{dd.MM.yyyy HH:mm}D %f %f %f %f %f %f %f',s_d,e_d);
s.magnetic.date = dataArray{1};
s.magnetic.gX = dataArray{3};
s.magnetic.gY = dataArray{4};
s.magnetic.gZ = dataArray{5};
s.magnetic.mX = dataArray{6};
s.magnetic.mY = dataArray{7};
s.magnetic.mZ = dataArray{8};

% % Align end date for pressure, acceleration and light
% max_date = max([ s.pressure.date(end) s.acceleration.date(end) s.light.date(end)]);
% 
% s.pressure.date = s.pressure.date(1):diff(s.pressure.date(1:2)):max_date;
% s.pressure.obs = [s.pressure.obs ; nan(numel(s.pressure.date)-numel(s.pressure.obs),1)];
% s.acceleration.date = s.acceleration.date(1):diff(s.acceleration.date(1:2)):max_date;
% s.acceleration.pit = [s.acceleration.pit ; nan(numel(s.acceleration.date)-numel(s.acceleration.pit),1)];
% s.acceleration.act = [s.acceleration.act ; nan(numel(s.acceleration.date)-numel(s.acceleration.act),1)];
% s.light.date = (s.light.date(1):diff(s.light.date(1:2)):max_date)';
% s.light.obs = [s.light.obs ; nan(numel(s.light.date)-numel(s.light.obs),1)];
end
 
function dataArray = readFile(filename,formatSpec,s_d,e_d,varargin)

if nargin>4
    plotit = varargin{1};
else
    plotit = false;
end

fileID = fopen(filename,'r');
if fileID==-1
    warning(['No file found for: ' filename])
    dataArray = cell(1,8);
    dataArray(1:8) = {NaN};
    dataArray(1) = {NaT};
    return
end
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'HeaderLines' ,6);
fclose(fileID);

% convert date to international standard
dataArray{1}.Format='default';

% Check if the date is on a regular grid
if numel(unique(diff(dataArray{1})))~=1
    % If not, then, it's a bit more tricky. 
    warning(['Date is not a regular grid for ' filename])
    % save old dataArray
    dataArrayOld = dataArray;
    % figure; histogram(diff(dataArray{1}))
    % Get median delta time
    dt = median(diff(dataArrayOld{1}));
    % figure; histogram(mod(datenum(dataArray{1}),datenum(dt)))
    % Define the new datetime as a regular grid with median date
    dataArray{1} = (s_d:dt:e_d)';
    for i=2:numel(dataArrayOld)
        dataArray{i}=interp1(dataArrayOld{1},dataArrayOld{i}, dataArray{1}, 'nearest','extrap');
        % figure; hold on; plot(dataArrayOld{1},dataArrayOld{i}); plot(dataArray{1},dataArray{i})
    end
else
    % if so, simply trim with stard and end date all values
    id_d = dataArray{1}>=s_d & dataArray{1}<e_d;
    
    if plotit
        figure; hold on; title(filename)
        plot(dataArray{1},dataArray{2})
        plot(dataArray{1}(id_d),dataArray{2}(id_d))
    end
    
    for i=1:numel(dataArray)
        dataArray{i}=dataArray{i}(id_d);
    end

end

end
 