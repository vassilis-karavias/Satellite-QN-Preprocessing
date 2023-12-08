% Simulation times

startTime = datetime(2023, 7,17,20,00,0);
stopTime = startTime + days(1);
sampleTime = 10;

% Satellite Scenario Initialisation
sc = satelliteScenario(startTime, stopTime, sampleTime);
satelliteScenarioViewer(sc);


% Satellite 1

semiMajorAxis = 6937890;
eccentricity = 0;
inclination = 97.7;
rightAscensionOfAscendingNode1 = 110;         % degrees
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
rightAscensionOfAscendingNode2 = 130;         % degrees
argumentOfPeriapsis = 0;                   % degrees
trueAnomaly = 0;                           % degrees

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
name = ["London", "Madrid", "Edinburgh", "Barcelona", "Paris", "Nantes", "Lyon", "Lisbon", "Bern", "Freiburg", "Milan", "Florence", "Naples", "Frankfurt", "Dortmund", "Berlin", "Athens", "Ioannina"];
lat = [51.5072, 40.4168, 55.9533, 41.3874, 48.8566, 47.2184, 45.7640, 38.7223, 46.9480, 47.9990, 45.4642, 43.7696, 40.8518, 50.1109, 51.5136, 52.5200, 37.9838, 39.6650];
lon = [-0.1276, -3.7038, -3.188267, 2.1686, 2.3522, -1.5538, 4.8357, -9.1393, 7.4474, 7.8421,9.19000, 11.2558, 14.2681, 8.6921, 7.4653, 13.4050, 23.7275, 20.8537];
alt = [11, 660, 47, 12, 35, 20, 237, 2, 540, 280, 120, 50, 100, 110, 86, 35, 20, 480];
gs = groundStation(sc, "Name", name, "Latitude",lat,"Longitude",lon, "Altitude",alt);


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
for groundstation = 1:size(gs,2)
    currentTime = [];
    currentTime2 = [];
    elevationAngle = [];
    elevationAngle2 = [];
    for t = 0:total_iterations
        time = startTime + seconds(t * sampleTime);
        
        [az,elev,r] = aer(gs(groundstation),sat1,time);
        if elev > 0
            currentTime(end + 1) = t * sampleTime;
            elevationAngle(end+1) = elev;
        end
        [az,elev2,r] = aer(gs(groundstation),sat2,time);
        if elev2 > 0
            currentTime2(end + 1) = t * sampleTime;
            elevationAngle2(end+1) = elev2;
        end
    end
    elevAngleTable = table;
    elevAngleTable.currentTime = transpose(currentTime);
    elevAngleTable.elevationAngle = transpose(elevationAngle);
    elevAngleTable2 = table;
    elevAngleTable2.currentTime = transpose(currentTime2);
    elevAngleTable2.elevationAngle = transpose(elevationAngle2);
    % elevAngleTable = table(transpose(currentTime), transpose(elevationAngle));
    elevAngleTable.currentDate = startTime + seconds(elevAngleTable.currentTime);
    elevAngleTable2.currentDate = startTime + seconds(elevAngleTable2.currentTime);
    writetable(elevAngleTable, strcat("satelliteElevAngle_", gs(groundstation).Name, "_",int2str(rightAscensionOfAscendingNode1),".csv"), "Delimiter",",");
    writetable(elevAngleTable2, strcat("satelliteElevAngle_", gs(groundstation).Name, "_",int2str(rightAscensionOfAscendingNode2),".csv"), "Delimiter",",");
end
% Play the Model 
play(sc);