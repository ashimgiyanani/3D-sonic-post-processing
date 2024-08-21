%% settings
scheme          = 'https';
host            = 'onedas.iwes.fraunhofer.de';
port         	= 443;
username        = 'test@onedas.org';
password        = '#test0/User1'; % password = input('Please enter your password: ')

dateTimeBegin 	= datetime(2020, 02, 01, 0, 0, 0, 'TimeZone', 'UTC');
dateTimeEnd 	= datetime(2020, 02, 02, 0, 0, 0, 'TimeZone', 'UTC');

% must all be of the same sample rate
channelPaths = { ...
    '/IN_MEMORY/TEST/ACCESSIBLE/T1/1 s_mean'
    '/IN_MEMORY/TEST/ACCESSIBLE/V1/1 s_mean'
};

%% load data
connector = OneDasConnector(scheme, host, port, username, password);
% without authentication: connector = OneDasConnector(scheme, host, port)

params.ChannelPaths = channelPaths;
data                = connector.Load(dateTimeBegin, dateTimeEnd, params);

%% prepare plot
y1          = data('/IN_MEMORY/TEST/ACCESSIBLE/T1/1 s_mean');
y2          = data('/IN_MEMORY/TEST/ACCESSIBLE/V1/1 s_mean');

sampleRate  = 1; % 1 Hz (adapt to your needs)
dt          = 1 / sampleRate / 86400;
time        = (datenum(dateTimeBegin) : dt : datenum(dateTimeEnd) - dt).';

%% plot
yyaxis left
plot(time, y1.Values)
ylabel([y1.Name ' / ' y1.Unit])
ylim([0 30])

yyaxis right
plot(time, y2.Values)
ylabel([y2.Name ' / ' y2.Unit])
ylim([0 30])

axis    = gca;
r1      = axis.YAxis(1);
r2      = axis.YAxis(2);
linkprop([r1 r2], 'Limits');

xlabel('Time')
xlim([time(1), time(end)])

title('OneDAS Connector Sample')
datetick('x', 'dd.mm HH:MM')
grid('on')