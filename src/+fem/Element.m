%% Element Class
%
%% Description
%
% This is an abstract super-class that generically specifies an
% <element.html element> shape in the FEMOOLab program.
% Essentially, this super-class declares abstract methods that define the
% general behavior of elements. These abstract methods are the functions
% that should be implemented in a derived sub-class that deals with
% specific types of elements.
%
%% Subclasses
%
% * <element_isoparametric.html Element_Isoparametric
% * <element_isogeometric.html Element_Isogeometric
%
%% Class definition
%
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General
        id    int32 = int32.empty;              % identification number
        gle   int32 = int32.empty;              % gather vector (stores element global d.o.f. numbers)
        anm   = [];                             % object of Anm class
        shape = [];                             % object of Shape class
        gauss = [];                             % object of Gauss class
        
        % Gauss integration quadrature
        gsystem_order int32  = int32.empty;     % order of Gauss quadrature for computation of global system arrays
        gderive_order int32  = int32.empty;     % order of Gauss quadrature for computation of derived variables
        gderive_npts  int32  = int32.empty;     % number of Gauss points for computation of derived variables        
        TGN           double = double.empty;    % transformation matrix of Gauss points results to nodal results
        
        % Physical attributes
        mat fem.Material = fem.Material.empty;  % object of Material class
        thk double       = double.empty;        % thickness
        
        % Domain forcing source
        src double = double.empty;              % components of body force (structural) / internal heat generation (thermal)
        
        % Natural boundary conditions (uniformly distributed over edges)
        lineNBC1 double = double.empty;         % matrix of standard NBCs (constant values) [corner1,corner2,forcing_components]
        lineNBC2 double = double.empty;         % matrix of radiative NBCs (d.o.f. dependent values) [corner1,corner2,forcing_term]
        
        % Convection conditions
        avgV   double = double.empty            % cartesian components of average nodal velocities
        normV  double = double.empty            % norm of average nodal velocities vector
        peclet double = double.empty            % Peclet number
        alpha  double = double.empty            % stabilization coefficient alpha
        beta   double = double.empty            % stabilization coefficient beta
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element()
            return;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Compute transformation matrix of gauss-to-node results.
        % Refs.:
        % -Hinton & Campbell, "Local and Global Smoothing of Discontinous
        % Finite Element Functions using a Least Squares Method",
        % Int. J. Num. Meth. Engng., Vol. 8, pp. 461-480, 1974.
        % -Burnett, D.S., "Finite Element Analysis - From Concepts to Applications",
        % Addison-Wesley, 1987.
        % -Martha, L.F., "Notas de Aula do Curso CIV 2118 - Metodo dos Elementos
        % Finitos", 1994.
        TGNmtx(this);

        %------------------------------------------------------------------
        % Initialize convection properties: velocities, Peclet number, and
        % stabilization coefficients.
        convProps(this);
        
        %------------------------------------------------------------------
        % Compute gradient matrix in a position of element domain given by
        % parametric coordinates.
        B = BmtxElem(this,J,r,s)
        
        %------------------------------------------------------------------
        % Compute Jacobian matrix in a position of element domain given by
        % parametric coordinates.
        J = JmtxDomain(this,r,s)
        
        %------------------------------------------------------------------
        % Compute Jacobian matrix in a position of element edge given by
        % parametric coordinates.
        J = JmtxEdge(this,edgLocIds,r)
        
        %------------------------------------------------------------------
        % Compute number of edge nodes and assemble vector of edge nodes ids
        [nedgen,edgLocIds] = edgeIds(this,corner1,corner2)
        
        %------------------------------------------------------------------
        % Assemble edge gather vector (stores local d.o.f.'s numbers).
        gledge = gleEdgeVct(this,nedgen,edgLocIds)
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % diffusive term: [B]'[C][B]h
        K = stiffDiffMtx(this)
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % convective term: [N]'[v]'[B]
        K = stiffConvMtx(this)
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % stabilization of convective term: b[B]'[V][B]
        % (steady-state analysis by SUPG method)
        K = stiffStabMtx(this)
        
        %------------------------------------------------------------------
        % Numerical integration of stiffness matrix part accounting for the
        % radiative boundary conditions over edges: [N]'[N]h
        K = stiffRadMtx(this)
        
        %------------------------------------------------------------------
        % Numerical integration of mass matrix: [N]'[N]m
        M = massMtx(this)
        
        %------------------------------------------------------------------
        % Numerical integration of equivalent nodal forcing vector from
        % standard natural boundary conditions prescribed over edges: [N]'q
        F = edgeEquivForceVct(this,~)
        
        %------------------------------------------------------------------
        % Numerical integration of equivalent nodal forcing vector from
        % internal domain source: [N]'Q
        F = domainEquivForceVct(this)
        
        %------------------------------------------------------------------
        % Numerical integration of nodal forcing vector accounting for the
        % stabilization of convective term: b[B]'[v]'Q
        % (steady-state analysis by SUPG method)
        F = domainStabForceVct(this)
        
        %------------------------------------------------------------------
        % Compute derived variables at Gauss points and gauss point cartesian
        % coordinates for a given element.
        [ngp,dvar,gpc] = derivedVar(this,u)
    end
end