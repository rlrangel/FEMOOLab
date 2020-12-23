classdef Print_Axisymmetric < Print
    %% Public exclusive methods
    % Exclusive methods of the Print_Axisymmetric sub-class.
    methods
        %------------------------------------------------------------------
        % Prints header of analysis results.
        function print = header(print,fid)
            fprintf(fid, '\n=================================================================\n');
            fprintf(fid, ' FEMOOLab - Finite Element Method Object-Oriented Laboratory\n');
            fprintf(fid, '    PONTIFICAL CATHOLIC UNIVERSITY OF RIO DE JANEIRO\n');
            fprintf(fid, '    DEPARTMENT OF CIVIL AND ENVIRONMENTAL ENGINEERING\n');
            fprintf(fid, '                          AND\n');
            fprintf(fid, '               TECGRAF/PUC-RIO INSTITUTE\n');
            fprintf(fid, '=================================================================\n');
        end
        
        %------------------------------------------------------------------
        % Prints analyis model type.
        function print = analysisLabel(print,fid)
            fprintf(fid, '\n\n----------------------------\n');
            fprintf(fid, 'ANALYSIS MODEL: ');
            fprintf(fid, 'AXISYMMETRIC\n');
            fprintf(fid, '----------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints global model description.
        function print = modelDescrip(print,model,fid)
            fprintf(fid, '\n\n----------------------------------------\n' );
            fprintf(fid, 'M O D E L  D E S C R I P T I O N\n' );
            fprintf(fid, '----------------------------------------\n');
            fprintf(fid, 'NUMBER OF NODES....................:%4d\n', model.nnp);
            fprintf(fid, 'NUMBER OF ELEMENTS ................:%4d\n', model.nel);
            fprintf(fid, 'NUMBER OF DEGREES OF FREEDOM.......:%4d\n', model.neq);
            fprintf(fid, 'NUMBER OF FREE DEGREES OF FREEDOM..:%4d\n', model.neqfree);
            fprintf(fid, 'NUMBER OF FIXED DEGREES OF FREEDOM.:%4d\n', model.neqfixed);
            fprintf(fid, 'NUMBER OF MATERIALS................:%4d\n', model.nmat);
            fprintf(fid, '----------------------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints material properties.
        function print = material(print,model,fid)
            fprintf(fid, '\n\n-------------------------------------\n');
            fprintf(fid, 'M A T E R I A L  P R O P E R T I E S\n');
            fprintf(fid, '-------------------------------------\n');
            if model.nmat > 0
                fprintf(fid, ' MATERIAL        E            G        POISSON\n');
                
                for m = 1:model.nmat
                    fprintf(fid, '%4d   %13.0f   %10.0f   %8.2f\n', m,...
                            model.materials(m).elasticity,...
                            model.materials(m).shear,...
                            model.materials(m).poisson);
                end
            else
                fprintf(fid, ' NO MATERIAL\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal coordinates.
        function print = nodalCoords(print,model,fid)
            fprintf(fid, '\n\n--------------------------------\n');
            fprintf(fid, 'N O D A L  C O O R D I N A T E S \n');
            fprintf(fid, '--------------------------------\n');
            fprintf(fid, ' NODE       COORD X          COORD Y\n');
            
            for n = 1:model.nnp
                fprintf(fid, '%4d   %11.3f   %14.3f\n', n,...
                        model.nodes(n).coords(1),...
                        model.nodes(n).coords(2));
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal restraint conditions.
        function print = nodalSupport(print,model,fid)
            fprintf(fid, '\n\n------------------------------\n');
            fprintf(fid, 'N O D A L  R E S T R A I N T S \n');
            fprintf(fid, '------------------------------\n');
            fprintf(fid, ' NODE   DISPL X   DISPL Y\n');
            
            for n = 1:model.nnp
                if(model.nodes(n).ebc(1) == 1)
                    node_restr1 = 'FIXED';
                else
                    node_restr1 = 'FREE';
                end
                
                if(model.nodes(n).ebc(2) == 1)
                    node_restr2 = 'FIXED';
                else
                    node_restr2 = 'FREE';
                end
                
                fprintf(fid, '%4d     %5s     %5s\n', n, node_restr1, node_restr2);
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal prescribed displacements.
        function print = nodalPrescDisp(print,model,fid)
            fprintf(fid, '\n\n--------------------------------\n');
            fprintf(fid, 'N O D A L  P R E S C.  D I S P L. \n');
            fprintf(fid, '--------------------------------\n');

            n_nprescdispl = 0;
            for n = 1:model.nnp
                if sum(model.nodes(n).prescDispl ~= 0) > 0
                    n_nprescdispl = n_nprescdispl + 1;
                end
            end

            if n_nprescdispl ~= 0
                fprintf(fid, ' NODE     DX          DY\n');

                for n = 1:model.nnp
                    if sum(model.nodes(n).prescDispl ~= 0) > 0
                        fprintf(fid, '%4d   %8.1f   %14.1f\n', n,...
                                model.nodes(n).prescDispl(1),...
                                model.nodes(n).prescDispl(2));
                    end
                end
            else
                fprintf(fid, ' NO PRESCRIBED DISPLACEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal loads.
        function print = nodalLoads(print,model,fid)
            fprintf(fid, '\n\n--------------------\n');
            fprintf(fid, 'N O D A L  L O A D S \n');
            fprintf(fid, '--------------------\n');

            n_nodalload = 0;
            for n = 1:model.nnp
                if sum(model.nodes(n).nodalLoad ~= 0) > 0
                    n_nodalload = n_nodalload + 1;
                end
            end

            if n_nodalload ~= 0
                fprintf(fid, ' NODE     FX           FY\n');

                for n = 1:model.nnp
                    if sum(model.nodes(n).nodalLoad ~= 0) > 0
                        fprintf(fid, '%4d    %9.3f    %12.3f\n', n,...
                                model.nodes(n).nodalLoad(1),...
                                model.nodes(n).nodalLoad(2));
                    end
                end
            else
                fprintf(fid, ' NO NODAL LOAD\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints elements information.
        function print = elements(print,model,fid)
            include_constants;

            fprintf(fid, '\n\n---------------\n');
            fprintf(fid, 'E L E M E N T S\n');
            fprintf(fid, '---------------\n');

            if model.nel > 0
                fprintf(fid, ' ELEMENT       TYPE       MAT     NODES\n');
            
                for e = 1:model.nel
                    mat = model.elems(e).material.id;
                    
                    element_type = model.elems(e).elem_type;
                    switch element_type
                        case TRI3
                            type = 'T3 element';
                            
                            n1 = model.elems(e).nodes(1).id;
                            n2 = model.elems(e).nodes(2).id;
                            n3 = model.elems(e).nodes(3).id;
                    
                            fprintf(fid, '%5d   %15s % 4d  %5d %d %d\n', ...
                                    e,type,mat,n1,n2,n3);
                                    
                        case QUAD4
                            type = 'Q4 element';
                            
                            n1 = model.elems(e).nodes(1).id;
                            n2 = model.elems(e).nodes(2).id;
                            n3 = model.elems(e).nodes(3).id;
                            n4 = model.elems(e).nodes(4).id;
                    
                            fprintf(fid, '%5d   %15s % 4d  %5d %d %d %d\n', ...
                                    e,type,mat,n1,n2,n3,n4);
                                
                        case TRI6
                            type = 'T6 element';
                            
                            n1 = model.elems(e).nodes(1).id;
                            n2 = model.elems(e).nodes(2).id;
                            n3 = model.elems(e).nodes(3).id;
                            n4 = model.elems(e).nodes(4).id;
                            n5 = model.elems(e).nodes(5).id;
                            n6 = model.elems(e).nodes(6).id;
                    
                            fprintf(fid, '%5d   %15s % 4d  %5d %d %d %d %d %d\n', ...
                                    e,type,mat,n1,n2,n3,n4,n5,n6);
                                
                        case QUAD8
                            type = 'Q8 element';
                            
                            n1 = model.elems(e).nodes(1).id;
                            n2 = model.elems(e).nodes(2).id;
                            n3 = model.elems(e).nodes(3).id;
                            n4 = model.elems(e).nodes(4).id;
                            n5 = model.elems(e).nodes(5).id;
                            n6 = model.elems(e).nodes(6).id;
                            n7 = model.elems(e).nodes(7).id;
                            n8 = model.elems(e).nodes(8).id;
                    
                            fprintf(fid, '%5d   %15s % 4d  %5d %d %d %d %d %d %d %d\n', ...
                                    e,type,mat,n1,n2,n3,n4,n5,n6,n7,n8);
                    end
                end
            else
                fprintf(fid, ' NO ELEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        function print = edgeLoads(print,model,fid)
            fprintf(fid, '\n\n----------------------------\n');
            fprintf(fid, 'E D G E  L O A D S\n');
            fprintf(fid, '----------------------------\n');

            n_load  = 0;
            for e = 1:model.nel
                if sum(model.elems(e).load.edgeLoad ~= 0) > 0
                    n_load = n_load + 1;
                end
            end

            if n_load ~= 0
                fprintf(fid, ' ELEMENT   NODE1   NODE2      QX         QY\n');
                
                for e = 1:model.nel
                    if sum(model.elems(e).load.edgeLoad ~= 0) > 0
                        n1 = model.elems(e).load.edgeLoad(1);
                        n2 = model.elems(e).load.edgeLoad(2);
                        qx = model.elems(e).load.edgeLoad(3);
                        qy = model.elems(e).load.edgeLoad(4);

                        fprintf(fid, '%5d   %6d %7d %11.3f %11.3f\n',...
                                e, n1, n2, qx, qy);
                    end
                end
            else
                fprintf(fid, ' NO EDGE LOAD\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        function print = areaLoads(print,model,fid)
            fprintf(fid, '\n\n----------------------------\n');
            fprintf(fid, 'A R E A  L O A D S\n');
            fprintf(fid, '----------------------------\n');

            n_load  = 0;
            for e = 1:model.nel
                if sum(model.elems(e).load.areaLoad ~= 0) > 0
                    n_load = n_load + 1;
                end
            end

            if n_load ~= 0
                fprintf(fid, ' ELEMENT       QX         QY\n');
                
                for e = 1:model.nel
                    if sum(model.elems(e).load.edgeLoad ~= 0) > 0
                        qx = model.elems(e).load.edgeLoad(1);
                        qy = model.elems(e).load.edgeLoad(2);

                        fprintf(fid, '%5d  %8.3f %10.3f\n',...
                                e, qx, qy);
                    end
                end
            else
                fprintf(fid, ' NO AREA LOAD\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of nodal displacement/rotation.
        function print = nodalDisplRot(print,model,fid)
            fprintf(fid, '\n\n-------------------------------------------\n');
            fprintf(fid, 'N O D A L  D I S P L A C E M E N T S\n');
            fprintf(fid, '------------------------------------------\n');
            fprintf(fid, ' NODE       DISPL X         DISPL Y\n');

            for n = 1:model.nnp
                dx = model.D(model.ID(1,n));
                dy = model.D(model.ID(2,n));
                fprintf(fid, '%4d     %11.6f     %11.6f\n', n, dx, dy);
            end
        end

        %------------------------------------------------------------------
        % Prints results of support reactions.
        function print = reactions(print,model,fid)
            fprintf(fid, '\n\n---------------------------------\n');
            fprintf(fid, 'S U P P O R T  R E A C T I O N S\n');
            fprintf(fid, '---------------------------------\n');
            fprintf(fid, ' NODE       FORCE X       FORCE Y\n');

            for n = 1:model.nnp
                if (model.ID(1,n) > model.neqfree) || (model.ID(2,n) > model.neqfree)

                    if model.ID(1,n) > model.neqfree 
                        reaction1 = model.F(model.ID(1,n));
                    else
                        reaction1 = 0.0;
                    end

                    if model.ID(2,n) > model.neqfree 
                        reaction2 = model.F(model.ID(2,n));
                    else
                        reaction2 = 0.0;
                    end

                    fprintf(fid, '%4d     %9.3f   %12.3f\n', n,reaction1,reaction2);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of stresses at gauss points.
        function print = gaussStresses(print,model,fid)
            fprintf(fid, '\n\n------------------------------\n');
            fprintf(fid, 'G A U S S  S T R E S S E S\n');
            fprintf(fid, '------------------------------\n');
            fprintf(fid, ' SIGMA_X       SIGMA_Y       TAU_XY\n');

            for e = 1:model.nel
                fprintf(fid, 'Elem %d\n',e);
                for n = 1:model.elems(e).n_gaussstress_pts
                    sx = model.sx_gp(n,e);
                    sy = model.sy_gp(n,e);
                    txy = model.txy_gp(n,e);
                    fprintf(fid, '%7.3f %13.3f %13.3f\n',sx,sy,txy);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of extrapolated nodal stress components.
        function print = elemExtrapStresses(print,model,fid)
            fprintf(fid, '\n\n------------------------------------\n');
            fprintf(fid, 'E X T R A P.  N O D A L  S T R E S S\n');
            fprintf(fid, '------------------------------------\n');
            fprintf(fid, ' SIGMA_X       SIGMA_Y       TAU_XY\n');

            for e = 1:model.nel
                fprintf(fid, 'Elem %d\n',e);
                for n = 1:model.elems(e).nen
                    sx = model.sx_elemextrap(n,e);
                    sy = model.sy_elemextrap(n,e);
                    txy = model.txy_elemextrap(n,e);
                    fprintf(fid, '%7.3f %13.3f %13.3f\n',sx,sy,txy);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of extrapolated node smooth stress components.
        function print = nodeExtrapStresses(print,model,fid)
            fprintf(fid, '\n\n----------------------------------\n');
            fprintf(fid, 'S M O O T H  N O D A L  S T R E S S\n');
            fprintf(fid, '----------------------------------\n');
            fprintf(fid, ' NODE     SIGMA_X       SIGMA_Y       TAU_XY\n');

            for n = 1:model.nnp
                sx =  model.sx_nodeextrap(n);
                sy =  model.sy_nodeextrap(n);
                txy = model.txy_nodeextrap(n);
                fprintf(fid, '%5d %10.3f %13.3f %13.3f\n',n,sx,sy,txy);
            end
        end
    end
    
    %% Public abstract methods
    % Implementation of the abstract methods declared in super-class <print.html *Print*>.
    methods
        %------------------------------------------------------------------
        % Prints analyis results.
        function results(print,model,fid)
            print.header(fid);
            fprintf(fid, '\n\n\n____________ M O D E L  I N F O R M A T I O N ____________\n');
            print.analysisLabel(fid);
            print.modelDescrip(model,fid);
            print.material(model,fid);
            print.nodalCoords(model,fid);
            print.nodalSupport(model,fid);
            print.nodalPrescDisp(model,fid);
            print.nodalLoads(model,fid);
            print.elements(model,fid);
            print.edgeLoads(model,fid);
            print.areaLoads(model,fid);
            fprintf(fid, '\n\n\n\n_____________ A N A L Y S I S  R E S U L T S _____________\n');
            print.nodalDisplRot(model,fid);
            print.reactions(model,fid);
            print.gaussStresses(model,fid);
            print.elemExtrapStresses(model,fid);
            print.nodeExtrapStresses(model,fid);
        end
    end
end