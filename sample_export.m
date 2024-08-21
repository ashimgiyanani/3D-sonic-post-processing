%% settings
scheme          = 'https';
host            = 'onedas.iwes.fraunhofer.de';
port         	= 443;
username        = 'test@onedas.org';
password        = '#test0/User1'; % password = input('Please enter your password: ')
targetFolder    = 'data';

dateTimeBegin 	= datetime(2020, 02, 01, 0, 0, 0, 'TimeZone', 'UTC');
dateTimeEnd 	= datetime(2020, 02, 02, 0, 0, 0, 'TimeZone', 'UTC');

% must all be of the same sample rate
channelPaths = { ...
    '/IN_MEMORY/TEST/ACCESSIBLE/T1/1 s_mean'
    '/IN_MEMORY/TEST/ACCESSIBLE/V1/1 s_mean'
};

%% export data
connector = OneDasConnector(scheme, host, port, username, password);
% without authentication: connector = OneDasConnector(scheme, host, port)

params.FileGranularity  = 'Day';    % Minute_1, Minute_10, Hour, Day, SingleFile
params.FileFormat       = 'MAT73';  % CSV, FAMOS, MAT73
params.ChannelPaths     = channelPaths;

connector.Export(dateTimeBegin, dateTimeEnd, params, targetFolder);