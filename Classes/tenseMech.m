classdef tenseMech<TensegritySettings    
    %Veřejné metody
    properties(Access = public)
        Property1
    end
    %Privátní metody
    properties(Access = private)
        
    end
    %Konstanty
    properties(Access = public,Constant)
        
    end
    
    methods
        function obj = tenseMech()
            obj@TensegritySettings()
        end
    end
end

