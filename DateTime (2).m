% DateTime - Used to format datetimes for display

% Copyright 2016 The MathWorks, Inc.

classdef (Hidden) DateTime < parallel.internal.display.DisplayableItem
    properties (SetAccess = immutable, GetAccess = private)
        Value
    end
    
    methods
        function obj = DateTime(displayHelper, value)
            obj@parallel.internal.display.DisplayableItem(displayHelper);
            obj.Value = value;
        end
        
        function displayInMATLAB(obj, name)
            obj.DisplayHelper.displayProperty(name, obj.DisplayHelper.formatDateTime(obj.Value));
        end
    end
    
end
