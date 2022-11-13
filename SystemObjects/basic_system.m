classdef basic_system < matlab.System
    % untitled4 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        m = 1;
        k = 10
        b = 1
    end

    properties (DiscreteState)

    end

    % Pre-computed constants
    properties (Access = private)

    end

    methods (Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            R = TensegritySettings();
        end

        function ydd = stepImpl(obj,u,y,yd)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            ydd = -(obj.k*y+obj.b*yd+u)/obj.m;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end
