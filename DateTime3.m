function d = datetime(varargin)
%DATETIME Create an array of codistributed datetimes.
%   
%   D = DATETIME(DS,'InputFormat',INFMT)
%   D = DATETIME(DS,'InputFormat',INFMT,'Locale',LOCALE)
%   D = DATETIME(DS,'InputFormat',INFMT,'PivotYear',PIVOT,...)
%   D = DATETIME(DV)
%   D = DATETIME(Y,MO,D,H,MI,S)
%   D = DATETIME(Y,MO,D)
%   D = DATETIME(Y,MO,D,H,MI,S,MS)
%   D = DATETIME(X,'ConvertFrom',TYPE)
%   D = DATETIME(X,'ConvertFrom','epochtime','Epoch',EPOCH)
%   D = DATETIME(...,'Format',FMT)
%   D = DATETIME(...,'TimeZone',TZ,...)
%    
%   LIMITATION:
%   When creating DATETIME from the strings in the cell array DS,
%   always specify the input format INFMT, for correctness.
%    
%   See also DATETIME.
%   


%   Copyright 2016 The MathWorks, Inc.

d = codistributed.pDatetimeDurationCommon(@datetime, [3 6], varargin{:});
