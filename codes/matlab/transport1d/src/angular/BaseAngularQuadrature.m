%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Base Angular Quadrature
%
%   Author:         Michael W. Hackemack
%   Institution:    
%   Year:           2019
%
%   Description:    MATLAB base class to 
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BaseAngularQuadrature < handle
    properties (Access = public)
        NumberDirections
        Directions
        Weights
        OppInds
        NumberMoments
        PnOrder
        Harmonics
        M2D
        D2M
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                            Constructor Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function obj = BaseAngularQuadrature(varargin)
            % do nothing at this time
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function compute_harmonics(obj)
            % Allocate arrays
            obj.Harmonics = zeros(obj.NumberMoments,obj.NumberDirections);
            obj.M2D       = zeros(obj.NumberMoments,obj.NumberDirections);
            obj.D2M       = zeros(obj.NumberMoments,obj.NumberDirections);
            % Compute harmonics
            for i=0:obj.PnOrder
                leg = legendre(i,obj.Directions);
                for q=1:obj.NumberDirections
                    obj.Harmonics(i+1,q) = leg(1,q);
                    obj.D2M(i+1,q) = obj.Weights(q)*obj.Harmonics(i+1,q);
                    obj.M2D(i+1,q) = (2*i+1)*obj.Harmonics(i+1,q);
                end
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
end