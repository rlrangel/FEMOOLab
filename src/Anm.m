%% Anm (Analysis Model) Class
%
% This is an abstract super-class that generically specifies an analysis 
% model in the StAnOOP program.
%
% Essentially, this super-class declares abstract methods that define
% the particular behavior of an analysis model. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis models.
%
classdef Anm < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = 0;  % flag for type of analysis model
        ndof = 0;  % number of degrees-of-freedom per node
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm(type,ndof)
            anm.type = type;
            anm.ndof = ndof;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Assemble material constitutive matrix of a given element.
        C = Cmtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble strain-displacement matrix at a given position of
        % an element.
        B = Bmtx(anm,elem,GradNcar,r,s);
        
        %------------------------------------------------------------------
        % Compute stress components (sx, sy, txy) at a given point of an element.
        str = pointStress(anm,C,B,d);
    end
end