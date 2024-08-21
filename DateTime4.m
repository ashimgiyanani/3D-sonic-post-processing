function d = datetime(arg1, varargin)
%DATETIME create tall array of datetime from tall arrays
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
%   Limitations:
%   When creating DATETIME from the strings in the cell array DS,
%   always specify the input format INFMT, for correctness.
%
%   See also DATETIME.
        
%   Copyright 2015-2016 The MathWorks, Inc.

% At least validate the first argument
arg1 = tall.validateType(arg1, mfilename, {'numeric', 'cellstr', 'string'}, 1);
d = slicefun(@(varargin) datetime(varargin{:}), arg1, varargin{:});
d = setKnownType(d, 'datetime');
end
