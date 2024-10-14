%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

clc, clear, close all, beep off, format long g
set(0,'defaultTextInterpreter','latex')
warning('off','MATLAB:namelengthmaxexceeded');

%#ok<*MINV>

%% Select Region

REGION = questdlg('Select Region','Region','Niger','Ganges-Brahmaputra','Niger');
if strcmp(REGION,'Niger')
    % Niger ------------------
    REGION_LAT = [2,20];
    REGION_LON = [-15,15];
    % ------------------------
elseif strcmp(REGION,'Ganges-Brahmaputra')
    % ganges-brahmaputra -----
    REGION_LAT = [20,33];
    REGION_LON = [72,100];
    % ------------------------
end

%% Load Altimetric Data

data.flag = string();
files.main = string();
k = 1;

% DAHITI
[file_d,path_d,index_d] = uigetfile('*.nc*','Select Dahiti Virtual Station Data Files','MultiSelect','on');
if index_d == 1
    if iscell(file_d)
        for i = 1:size(file_d,2)
            ncid = netcdf.open([path_d,char(file_d(i))], 'NC_NOWRITE');
            files.main(k,1) = 'WL' + extractBefore(string(file_d(i)),'.nc');
            data.main.(files.main(k,1)).datetime = datetime(netcdf.getVar(ncid,0));
            data.main.(files.main(k,1)).water_level = netcdf.getVar(ncid,1);
            data.main.(files.main(k,1)).error = netcdf.getVar(ncid,2);
            data.main.(files.main(k,1)).latitude = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'latitude');
            data.main.(files.main(k,1)).longitude = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'longitude');
            data.flag(k,1) = 'Dahiti';
            k = k + 1;
        end
    else
        files.main(k,1) = extractBefore(string(file_d),'.nc');
        data.main.(files.main) = rncdf([path_wse,char(file_wse)]);
        data.main.(files.main).datetime = datetime(data.main.(files.main).datetime);
        data.flag(k,1) = 'Dahiti';
        k = k + 1;
    end % DAHITI
end

% HYDROWEB
[file_h,path_h,index_h] = uigetfile('*.txt*','Select Hydroweb Virtual Station Data Files','MultiSelect','on');
if index_h == 1
    if iscell(file_h)
        for i = 1:size(file_h,2)
            files.main(k,1) = extractBefore(string(file_h(i)),'.txt');
            if contains(files.main(k,1),'-')
                files.main(k,1) = strrep(files.main(k,1),'-','_');
            end
            data.main.(files.main(k,1)) = HWread([path_h,char(file_h(i))]);
            data.flag(k,1) = 'Hydroweb';
            k = k + 1;
        end
    else
        files.main(k,1) = extractBefore(string(file_h),'.nc');
        data.main.(files.main) = HWread([path_wse,char(file_wse)]);
        data.flag(k,1) = 'Hydroweb';
        k = k + 1;
    end % HYDROWEB
end

% CLMS
[file_c,path_c,index_c] = uigetfile('*.json*','Select CLMS Virtual Station Data Files','MultiSelect','on');
if index_c == 1
    if iscell(file_c)
        for i = 1:size(file_c,2)
            files.main(k,1) = extractBetween(string(file_c(i)),'c_gls_','_0');
            data.main.(files.main(k,1)) = jsondecode(fileread([path_c,char(file_c(i))]));
            data.main.(files.main(k,1)).identifier = string(vertcat(data.main.(files.main(k,1)).data.identifier));
            dt = vertcat(data.main.(files.main(k,1)).data.datetime); 
            dt = datetime([double(string(dt(:,1:4))), double(string(dt(:,6:7))), double(string(dt(:,9:10)))]);
            data.main.(files.main(k,1)).datetime = dt;
            data.main.(files.main(k,1)).wse = vertcat(data.main.(files.main(k,1)).data.orthometric_height_of_water_surface_at_reference_position);
            data.main.(files.main(k,1)).uncertainty = vertcat(data.main.(files.main(k,1)).data.associated_uncertainty);
            % data.main.(files.main(k,1)).satellite = string(vertcat(data.main.(files.main(k,1)).data.satellite));
            % data.main.(files.main(k,1)).ground_track_number = vertcat(data.main.(files.main(k,1)).data.ground_track_number);
            data.main.(files.main(k,1)) = rmfield(data.main.(files.main(k,1)),'data');
            data.flag(k,1) = 'CLMS';
            k = k + 1;
        end
    else
        files.main(k,1) = extractBetween(string(file_c),'c_gls_','_0');
        data.main.(files.main) = jsondecode(fileread([path_c,char(file_c)]));
        data.main.(files.main).identifier = string(vertcat(data.main.(files.main).data.identifier));
        dt = vertcat(data.main.(files.main).data.datetime); 
        dt = datetime([double(string(dt(:,1:4))), double(string(dt(:,6:7))), double(string(dt(:,9:10)))]);
        data.main.(files.main).datetime = dt;
        data.main.(files.main).wse = vertcat(data.main.(files.main).data.orthometric_height_of_water_surface_at_reference_position);
        data.main.(files.main).uncertainty = vertcat(data.main.(files.main).data.associated_uncertainty);
        % data.main.(files.main).satellite = string(vertcat(data.main.(files.main).data.satellite));
        % data.main.(files.main).ground_track_number = vertcat(data.main.(files.main).data.ground_track_number);
        data.main.(files.main) = rmfield(data.main.(files.main),'data');
        data.flag(k,1) = 'CLMS';
        k = k + 1;
    end % CLMS
end

clear dt

SOURCES = unique(data.flag);
if isempty(SOURCES)
    error('No Data From Any Resources Were Selected!')
end

%% Time series range

years.start = zeros(length(files.main),1);
years.end = zeros(length(files.main),1);
for i = 1:length(files.main)
    data_i =  data.main.(files.main(i));
    if strcmp(data.flag(i), 'Dahiti')
        years.start(i) = data_i.datetime(1).Year;
        years.end(i) = data_i.datetime(end).Year;
    elseif strcmp(data.flag(i), 'Hydroweb')
        years.start(i) = data_i.FirstDate.Year;
        years.end(i) = data_i.LastDate.Year;
    elseif strcmp(data.flag(i), 'CLMS')
        years.start(i) = data_i.datetime(1).Year;
        years.end(i) = data_i.datetime(end).Year;
    end
end

%% Sort Stations

[years.start,idx] = sort(years.start,'descend');
years.end = years.end(idx);
files.main = files.main(idx);
data.flag = data.flag(idx);
idx_h = data.flag == 'Hydroweb';
idx_d = data.flag == 'Dahiti';
idx_c = data.flag == 'CLMS';

%% Data Matrix

SY = min(years.start);
EY = max(years.end);

time1 = [repelem((SY:EY)',12), ...
        repmat((1:12)',length(SY:EY)',1), ...
        ones((EY-SY+1)*12,1)*15];
time1 = datetime(time1);
time = time1.Year + time1.Month/12 + time1.Day/365;
X1 = NaN(length(files.main),length(time1));

lat = zeros(length(files.main),1);
lon = zeros(length(files.main),1);
coords.vs.main = table(lat,lon);
mean_wl = zeros(length(files.main),1);
fprintf('\nCreating data matrix for visualizing time series...\n')
for i = 1:length(files.main)
    if strcmp(data.flag(i),'Dahiti')
        data_i = data.main.(files.main(i));
        for j = 1:length(time1)
            idx = (time1.Year(j) == data_i.datetime.Year) + (time1.Month(j) == data_i.datetime.Month);
            idx = idx == 2;
            X1(i,j) = mean(data_i.water_level(idx));
            coords.vs.main.lat(i) = data_i.latitude;
            coords.vs.main.lon(i) = data_i.longitude;
            mean_wl(i) = mean(data_i.water_level);
        end
    elseif strcmp(data.flag(i),'Hydroweb')
        data_i = data.main.(files.main(i)).MainData;
        for j = 1:length(time1)
            idx = (time1.Year(j) == data_i.Time.Year) + (time1.Month(j) == data_i.Time.Month);
            idx = idx == 2;
            X1(i,j) = mean(data_i.("ORTHOMETRIC HEIGHT (M) OF WATER SURFACE AT REFERENCE POSITION")(idx));
            tmp = data_i.("LATITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=90);
            coords.vs.main.lat(i) = mean(tmp);
            tmp = data_i.("LONGITUDE OF ALTIMETRY MEASUREMENT (deg)"); tmp = tmp(abs(tmp)<=180);
            coords.vs.main.lon(i) = mean(tmp);
            mean_wl(i) = mean(data_i.("ORTHOMETRIC HEIGHT (M) OF WATER SURFACE AT REFERENCE POSITION"));
        end
    elseif strcmp(data.flag(i),'CLMS')
        data_i = data.main.(files.main(i));
        for j = 1:length(time1)
            idx = (time1.Year(j) == data_i.datetime.Year) + (time1.Month(j) == data_i.datetime.Month);
            idx = idx == 2;
            X1(i,j) = mean(data_i.wse(idx));
            coords.vs.main.lat(i) = data_i.geometry.coordinates(2);
            coords.vs.main.lon(i) = data_i.geometry.coordinates(1);
            mean_wl(i) = mean(data_i.wse);
        end
    end
    if i == round(0.25*length(files.main))
        disp('Progress: 25%');
    elseif i == round(0.50*length(files.main))
        disp('Progress: 50%');
    elseif i == round(0.75*length(files.main))
        disp('Progress: 75%');
    elseif i == length(files.main)
        disp('Progress: 100%');
    end
end

%% Plots

figure()
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
ms = 100;
if any(ismember(SOURCES,'CLMS'))
    geoscatter(coords.vs.main.lat(idx_c),coords.vs.main.lon(idx_c),ms,mean_wl(idx_c),'filled','square')
end
if any(ismember(SOURCES,'Dahiti'))
    geoscatter(coords.vs.main.lat(idx_d),coords.vs.main.lon(idx_d),ms,mean_wl(idx_d),'filled','^')
end
if any(ismember(SOURCES,'Hydroweb'))
    geoscatter(coords.vs.main.lat(idx_h),coords.vs.main.lon(idx_h),ms,mean_wl(idx_h),'filled','o')
end

geolimits(REGION_LAT,REGION_LON)
clim([0, max(mean_wl)])
cb = colorbar; colormap('turbo')
set(get(cb,'ylabel'),'String','Mean Water Level [meters]','FontSize',12,'Interpreter','latex');
title(sprintf('%s - %g Virtual Stations',REGION,length(files.main)),'FontSize',13)
legend(SOURCES,'Interpreter','latex','FontSize',12)
gx = gca();
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';


figure()
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
ms = 15; hold on
if any(ismember(SOURCES,'CLMS'))
    geoplot(coords.vs.main.lat(idx_c),coords.vs.main.lon(idx_c),'or','MarkerSize',ms,'LineWidth',1.5)
end
if any(ismember(SOURCES,'Dahiti'))
    geoplot(coords.vs.main.lat(idx_d),coords.vs.main.lon(idx_d),'^b','MarkerSize',ms,'LineWidth',1.5)
end
if any(ismember(SOURCES,'Hydroweb'))
    geoplot(coords.vs.main.lat(idx_h),coords.vs.main.lon(idx_h),'xg','MarkerSize',ms,'LineWidth',1.5)
end

geolimits(REGION_LAT,REGION_LON)
title(sprintf('%s - %g Virtual Stations',REGION,length(files.main)),'FontSize',13)
legend(SOURCES,'Interpreter','latex','FontSize',12)
gx = gca();
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';


figure()
im = imagesc(time1.Year+time1.Month/12+time1.Day/365,(1:length(files.main))',X1);
clim([0 max(X1(:))])
title(sprintf('%s - Monthly Water Level Heatmap',REGION),'FontSize',13)
xlabel('Time','FontSize',12); ylabel('Station','FontSize',12)
[~,tmp] = unique(time1.Year); tmp = time1.Year(tmp);
xticks(tmp(1:2:end))
yticklabels('')
cmap = sky; cmap = [ones(1,3)*0.1; cmap]; colormap(cmap)
cb = colorbar; 
set(get(cb,'ylabel'),'String','Water Level [meters]','FontSize',12,'Interpreter','latex');
ax = gca();
ax.TickLabelInterpreter = 'latex';


%% ------------------------ Functions --------------------------------

function data = HWread(filepath)

data = rheader(filepath);
MainData = readtimetable(filepath,'NumHeaderLines',data.EOH,'TextType','string');
MainData.Var4 = [];
time = split(MainData.Var1,':');
dt = datetime(double([MainData.Time.Year,MainData.Time.Month,MainData.Time.Day,time(:,1),time(:,2),zeros(length(time),1)]));
MainData.Time = dt;
MainData.Var1 =[];
MainData.Properties.VariableNames = data.Content(3:end);
data.MainData = MainData;

function header = rheader(filepath)
    
        FileID = fopen(filepath);
        header.EOH = 0;
        while true

            header.EOH = header.EOH + 1;
            line = fgetl(FileID);

            if contains(line,'#############')
                fclose(FileID);
                break
            end

            if contains(line,'#BASIN')
                header.Basin = strtrim(line(9:end));
            end

            if contains(line,'#RIVER')
                header.River = strtrim(line(9:end));
            end

            if contains(line,'#ID')
                header.ID = strtrim(line(6:end));
            end

            if contains(line,'#TRIBUTARY OF')
                header.Tributary = strtrim(line(16:end));
            end

            if contains(line,'#APPROX. WIDTH OF REACH (m)')
                header.ApproxWidthOfReach = str2double(line(30:end));
            end

            if contains(line,'SURFACE OF UPSTREAM WATERSHED (km2)')
                header.SurfaceOfUpstreamWatershed = str2double(line(39:end));
            end

            if contains(line,'#RATING CURVE PARAMETERS')
                header.RatingCurveParameters = sscanf(line(67:end),'%f');
                if isempty(header.RatingCurveParameters)
                    header.RatingCurveParameters = NaN;
                end
            end

            if contains(line,'#REFERENCE ELLIPSOID')
                header.Refrence.Ellipsoid = strtrim(line(23:end));
            end

            if contains(line,'#REFERENCE LONGITUDE')
                header.Refrence.Longitude = str2double(line(23:end));
            end

            if contains(line,'#REFERENCE LATITUDE')
                header.Refrence.Latitude = str2double(line(22:end));
            end

            if contains(line,'#REFERENCE DISTANCE (km)')
                header.Refrence.Distance = str2double(line(27:end));
            end

            if contains(line,'#GEOID MODEL')
                header.GeoidModel = strtrim(line(15:end));
            end

            if contains(line,'#GEOID ONDULATION AT REF POSITION(M.mm)')
                header.GeoidOnulation = str2double(line(42:end));
            end

            if contains(line,'##MISSION(S)-TRACK(S)')
                header.MissionTrack = strtrim(line(23:end));
            end

            if contains(line,'#STATUS')
                header.Status = strtrim(line(10:end));
            end

            if contains(line,'#VALIDATION CRITERIA')
                header.ValidationCriteria = strtrim(line(23:end));
            end

            if contains(line,'#MEAN ALTITUDE(M.mm)')
                header.MeanAltitude = str2double(line(23:end));
            end

            if contains(line,'#MEAN SLOPE (mm/km)')
                header.MeanSlope = str2double(line(22:end));
            end

            if contains(line,'#NUMBER OF MEASUREMENTS IN DATASET')
                header.MeasurmentsNum = str2double(line(37:end));
            end

            if contains(line,'#FIRST DATE IN DATASET')
                header.FirstDate = datetime(sscanf(strrep(line(25:end),'-',' '),'%f')');
            end

            if contains(line,'#LAST DATE IN DATASET')
                header.LastDate = datetime(sscanf(strrep(line(24:end),'-',' '),'%f')');
            end

            if contains(line,'#DISTANCE MIN IN DATASET (km)')
                header.MinDistance = str2double(line(32:end));
            end

            if contains(line,'#DISTANCE MAX IN DATASET (km)')
                header.MaxDistance = str2double(line(32:end));
            end

            if contains(line,'#PRODUCTION DATE')
                header.ProductionDate = datetime(sscanf(strrep(line(19:end),'-',' '),'%f')');
            end

            if contains(line,'#PRODUCT VERSION')
                header.ProductionVersion = strtrim(line(19:end));
            end

            if contains(line,'#PRODUCT CITATION')
                header.ProductCitation = strtrim(line(20:end));
            end

            if contains(line,'#SOURCES')
                header.Sources = strtrim(line(11:end));
                if isempty(header.Sources)
                    header.Sources = NaN;
                end
            end

            if contains(line,'#COL')
                cnum = str2double(extractBefore(line(5:end),' : '));
                if cnum == 1
                    header.Content = "";
                end
                header.Content(cnum,1) = extractAfter(line,' : ');
            end

        end
    end

end