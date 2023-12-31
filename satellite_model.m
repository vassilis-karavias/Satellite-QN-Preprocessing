% Simulation times

startTime = datetime(2023, 7,17,20,00,0);
stopTime = startTime + days(1);
sampleTime = 1;

% Satellite Scenario Initialisation
sc = satelliteScenario(startTime, stopTime, sampleTime);
satelliteScenarioViewer(sc);


% Satellite 1

semiMajorAxis = 6937890;
eccentricity = 0;
inclination = 97.7;
rightAscensionOfAscendingNode1 = 90;         % degrees
argumentOfPeriapsis = 0;                   % degrees
trueAnomaly = 0;                           % degrees

sat1 = satellite(sc, ...
    semiMajorAxis, ...
    eccentricity, ...
    inclination, ...
    rightAscensionOfAscendingNode1, ...
    argumentOfPeriapsis, ...
    trueAnomaly, ...
    "Name","Satellite 1", ...
    "OrbitPropagator","two-body-keplerian");


% Satellite 2


semiMajorAxis = 6951000;
eccentricity = 0;
inclination = 97.5;
rightAscensionOfAscendingNode2 = 80;         % degrees
argumentOfPeriapsis = 0;                   % degrees
trueAnomaly = 0;                           % degrees
% 
sat2 = satellite(sc, ...
    semiMajorAxis, ...
    eccentricity, ...
    inclination, ...
    rightAscensionOfAscendingNode2, ...
    argumentOfPeriapsis, ...
    trueAnomaly, ...
    "Name","Satellite 2", ...
    "OrbitPropagator","two-body-keplerian");


% Add GroundTrack To Model


groundTrack(sat1, "LeadTime",12000);
groundTrack(sat2, "LeadTime",12000);


% Add Ground Stations
% name = ["London", "Madrid", "Edinburgh", "Barcelona", "Paris", "Nantes", "Lyon", "Lisbon", "Bern", "Freiburg", "Milan", "Florence", "Naples", "Frankfurt", "Dortmund", "Berlin", "Athens", "Ioannina"];
% lat = [51.5072, 40.4168, 55.9533, 41.3874, 48.8566, 47.2184, 45.7640, 38.7223, 46.9480, 47.9990, 45.4642, 43.7696, 40.8518, 50.1109, 51.5136, 52.5200, 37.9838, 39.6650];
% lon = [-0.1276, -3.7038, -3.188267, 2.1686, 2.3522, -1.5538, 4.8357, -9.1393, 7.4474, 7.8421,9.19000, 11.2558, 14.2681, 8.6921, 7.4653, 13.4050, 23.7275, 20.8537];
% alt = [11, 660, 47, 12, 35, 20, 237, 2, 540, 280, 120, 50, 100, 110, 86, 35, 20, 480];
% gs = groundStation(sc, "Name", name, "Latitude",lat,"Longitude",lon, "Altitude",alt);


name = ["London", "Graz", "Frankfurt", "Johanessburg", "Sao Paulo","Tokyo", "Auckland", "New Delhi", "Mumbai", "Bangalore", "Xinglong", "Nanshan", "Perth", "Brisbane", "Melborne", "Albany", "Denver"];
lat = [51.5072, 47.0707, 50.1109, -26.2041, -23.5558, 35.653, -36.8509, 28.6139, 19.0760, 12.9716, 40.395833, 43.475278, -31.9523, -27.4705, -37.8136, 42.6526, 39.7392];
lon = [-0.1276, 15.4395, 8.6821, 28.0473, -46.6496, 139.839, 174.7645, 77.2090, 72.8777, 77.5946, 117.5775, 87.176667, 115.8613, 153.0260, 144.9631, -73.7562, -104.9903];
alt = [11, 490, 112, 1753, 760, 40, 200, 216, 14, 920, 890, 2028, 0, 32, 31, 45, 1609];
gs = groundStation(sc, "Name", name, "Latitude",lat,"Longitude",lon, "Altitude",alt);
% 
% name = ["Nashville"];
% lat = [36.1627];
% lon = [-86.7816];
% alt = [169];
% gs = groundStation(sc, "Name", name, "Latitude",lat,"Longitude",lon, "Altitude",alt);
% get access times of satellite with gs:

% ac = access(sat1, gs);
% intvls = accessIntervals(ac);
% writetable(intvls, "satelliteoverpass_1.csv", "Delimiter",",")
% ac = access(sat2, gs);
% intvls = accessIntervals(ac);
% writetable(intvls, "satelliteoverpass_2.csv", "Delimiter",",")


% get azimuth angle, elevation andge and range at times:
total_time = datevec(stopTime - startTime);
total_time_sec = total_time(2) * 2592000 + total_time(3) * 86400 + total_time(4) * 3600 + total_time(5) * 60 + total_time(6);
total_iterations = floor(total_time_sec/sampleTime);
% 
for groundstation = 1:size(gs,2)
    currentTime = [];
    currentTime2 = [];
    elevationAngle = [];
    elevationAngle2 = [];
    range = [];
    range2 = [];
    for t = 0:total_iterations
        time = startTime + seconds(t * sampleTime);

        [az,elev,r] = aer(gs(groundstation),sat1,time);
        if elev > 0
            currentTime(end + 1) = t * sampleTime;
            elevationAngle(end+1) = elev;
            range(end+1) = r;
        end
        [az,elev2,r2] = aer(gs(groundstation),sat2,time);
        if elev2 > 0
            currentTime2(end + 1) = t * sampleTime;
            elevationAngle2(end+1) = elev2;
            range2(end+1) = r2;
        end
    end
    elevAngleTable = table;
    elevAngleTable.currentTime = transpose(currentTime);
    elevAngleTable.elevationAngle = transpose(elevationAngle);
    elevAngleTable.Range = transpose(range);
    elevAngleTable2 = table;
    elevAngleTable2.currentTime = transpose(currentTime2);
    elevAngleTable2.elevationAngle = transpose(elevationAngle2);
    elevAngleTable2.range = transpose(range2);
    % elevAngleTable = table(transpose(currentTime), transpose(elevationAngle));
    elevAngleTable.currentDate = startTime + seconds(elevAngleTable.currentTime);
    elevAngleTable2.currentDate = startTime + seconds(elevAngleTable2.currentTime);
    writetable(elevAngleTable, strcat("global_satelliteElevAngle_", gs(groundstation).Name, "_",int2str(rightAscensionOfAscendingNode1),".csv"), "Delimiter",",");
    writetable(elevAngleTable2, strcat("global_satelliteElevAngle_", gs(groundstation).Name, "_",int2str(rightAscensionOfAscendingNode2),".csv"), "Delimiter",",");
end
% Play the Model 
play(sc);