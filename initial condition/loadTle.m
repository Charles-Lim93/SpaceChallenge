%==========================================================================
% loadTle
%   load TLE data from given data file
%
%   Given:          
%
%   Returned (function value):
%     TLE                               structure       Two line element
%       satelliteNumber                 
%       classification
%       launchYear                      ( Year ) 
%       launchNumber                    
%       pieceOfLaunch
%       epochYear                       ( Year )
%       epochDay                        ( Day )
%       firstMeanMotionDerivative       ( 2 rev/day ) 
%       secondMeanMotionDerivative      ( 6 rev/day^2 ) 
%       dragTerm
%       ephemerisType   
%       elementNumber
%       checksum1
%       inclination                     ( deg )
%       ascendingNode                   ( deg )
%       eccentricity    
%       argumentOfPerigee               ( deg )
%       meanAnomaly                     ( deg )
%       meanMotion                      ( rev/day )
%       checksum2

%      
%   References:
%
%       http://celestrak.com/NORAD/documentation/tle-fmt.asp
%       (NORAD site)
%
%   This revision:  2011 March 25
%
%==========================================================================

function Tle = loadTle(path, tleFileName)

%     [tleFileName, tlePathName] = uigetfile('*.txt','Open TLE data'); 
    tlePathName = path;
    fid = fopen([tlePathName,'\',tleFileName]);
    
    % Zeroth line, satellite name
    Tle.satelliteName = fgetl(fid);
    
    % First line 
    tline = fgetl(fid);
    Tle.satelliteNumber = str2double( tline(3:7) );
    Tle.classification = tline(8);
    Tle.launchYear = str2double(tline(10:11));
    Tle.launchNumber = str2double(tline(12:14));
    Tle.pieceOfLaunch = tline(15:17);
    Tle.epochYear = str2double(tline(19:20));
    Tle.epochDay = str2double(tline(21:32));
    Tle.firstMeanMotionDerivative = str2double(tline(34:43));
    Tle.secondMeanMotionDerivative = ...
                   str2double(tline(45:50))*10^str2double(tline(51:52))...
                   *1e-5; 
    Tle.dragTerm = str2double(tline(54:59))*10^str2double(tline(60:61))...
                   *1e-5; 
    Tle.ephemerisType = tline(63);
    Tle.elementNumber = str2double(tline(65:68));
    Tle.checksum1 = str2double(tline(69));
        
    if Tle.launchYear < 57;  Tle.launchYear = Tle.launchYear + 2000; end;
    if Tle.epochYear < 57;  Tle.epochYear = Tle.epochYear + 2000; end;
    
    % Second line
    tline = fgetl(fid);
    Tle.inclination = str2double(tline(9:16));
    Tle.ascendingNode = str2double(tline(18:25));    
    Tle.eccentricity = str2double(tline(27:33))*1e-7;        
    Tle.argumentOfPerigee = str2double(tline(35:42));
    Tle.meanAnomaly = str2double(tline(44:51));    
    Tle.meanMotion = str2double(tline(53:63));
    Tle.checksum2 = str2double(tline(69)); 
   
    fclose(fid);  
    
end
