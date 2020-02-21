%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Mesh Class
%
%   Author:         Michael W. Hackemack
%   Institution:    
%   Year:           2019
%
%   Description:    MATLAB class to store 
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Mesh < handle
    properties (Access = public)
        Dimension = 1
        TotalVertices
        TotalCells
        TotalFaces
        TotalInteriorFaces
        TotalBoundaryFaces
        InteriorFaces
        BoundaryFaces
    end
    properties (Access = public)
        Vertices
        CellVerts
        CellFaces
        FaceCells
        FaceVerts
        MatID
        FaceID
        CellVertexNumbers
        VertexCells
        VertexFaces
        CellFaceVerts
        CellNeighbors
        CellNeighborFaces
    end
    properties (Access = public)
        CellCenter
        FaceCenter
        CellVolume
        CellSurfaceArea
        FaceArea
        FaceNormal
        OrthogonalProjection
    end
    properties (Access = public)
        MaxIrregularity = inf
        OriginalCellCount
        OriginalFaceCount
        OriginalVertexCount
        MeshRefinementLevel
        PreviousCell
        CellRefinedLastCycle
        CellRefinementFlag
        CellRefinementLevel
        CellRefinementTree
        CellRefinementTop
        CellRefinementTreeHierarchy
        CellRefinementTreeHierarchyLevel
    end
    properties (Access = public)
        x
        Lx
        minX
        maxX
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                            Constructor Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function obj = Mesh(varargin)
            if nargin == 1
                obj.x = varargin{1};
                obj.Constructor_Points();
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function Constructor_Points(obj)
            obj.Lx = length(obj.x);
            % Define Total Geometry Space
            if size(obj.x,1) < size(obj.x,2)
                obj.x = obj.x';
            end
            obj.Vertices = obj.x;
            obj.TotalVertices = length(obj.Vertices);
            obj.TotalCells = obj.TotalVertices - 1;
            obj.TotalFaces = obj.TotalVertices;
            obj.TotalInteriorFaces = obj.TotalFaces - 2;
            obj.TotalBoundaryFaces = 2;
            % Build Remaining Structures
            obj.Allocate_Arrays();
            for c=1:obj.TotalCells
                obj.CellVerts{c} = [c,c+1];
                obj.CellFaceVerts{c}{1} = c;
                obj.CellFaceVerts{c}{2} = c+1;
                obj.CellFaces{c} = [c,c+1];
                obj.CellVolume(c) = obj.Vertices(c+1) - obj.Vertices(c);
                obj.CellCenter(c) = (obj.Vertices(c+1) + obj.Vertices(c))/2;
                obj.FaceVerts{c} = c;
                if c==1
                    obj.CellNeighbors{c} = c+1;
                    obj.CellNeighborFaces{c} = c+1;
                    obj.CellSurfaceArea(c) = 1;
                elseif c==obj.TotalCells
                    obj.CellNeighbors{c} = c-1;
                    obj.CellNeighborFaces{c} = c;
                    obj.CellSurfaceArea(c) = 1;
                else
                    obj.CellNeighbors{c} = [c-1,c+1];
                    obj.CellNeighborFaces{c} = [c,c+1];
                    obj.CellSurfaceArea(c) = 2;
                end
            end
            obj.FaceVerts{end} = obj.TotalCells + 1;
            obj.FaceID = [1,zeros(1,obj.TotalFaces-2),1]';
            obj.FaceArea = ones(length(obj.FaceVerts),1);
            obj.FaceCenter = obj.Vertices;
            obj.FaceNormal = [-1,ones(1,obj.TotalFaces-1)]';
            for f=1:obj.TotalFaces
                if f==1
                    obj.OrthogonalProjection(1,1) = obj.x(2) - obj.x(1);
                    obj.VertexCells{1} = 1;
                    obj.FaceCells(1,1) = 1;
                elseif f==obj.TotalFaces
                    obj.OrthogonalProjection(end,1) = obj.x(end) - obj.x(end-1);
                    obj.VertexCells{obj.TotalFaces} = obj.Lx - 1;
                    obj.FaceCells(end,1) = obj.TotalCells;
                else
                    obj.OrthogonalProjection(f,1) = obj.x(f) - obj.x(f-1);
                    obj.OrthogonalProjection(f,2) = obj.x(f+1) - obj.x(f);
                    obj.VertexCells{f} = [f-1,f];
                    obj.FaceCells(f,1) = f-1;
                    obj.FaceCells(f,2) = f;
                end
            end
            obj.minX = min(obj.Vertices(:,1));
            obj.maxX = max(obj.Vertices(:,1));
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function Allocate_Arrays(obj)
            % Cell Arrays
            obj.MatID = ones(obj.TotalCells, 1);
            obj.CellVerts = cell(obj.TotalCells, 1);
            obj.CellFaceVerts = cell(obj.TotalCells, 1);
            obj.CellNeighbors = cell(obj.TotalCells, 1);
            obj.CellNeighborFaces = cell(obj.TotalCells, 1);
            obj.CellCenter = zeros(obj.TotalCells, obj.Dimension);
            obj.CellVolume = zeros(obj.TotalCells, 1);
            obj.CellSurfaceArea = zeros(obj.TotalCells, 1);
            obj.CellFaces = cell(obj.TotalCells, 1);
            % Face Arrays
            obj.FaceVerts = cell(obj.TotalFaces, 1);
            obj.FaceCells = zeros(obj.TotalFaces, 2);
            obj.InteriorFaces = zeros(obj.TotalInteriorFaces,1);
            obj.BoundaryFaces = zeros(obj.TotalBoundaryFaces,1);
            obj.OrthogonalProjection = zeros(obj.TotalFaces,2);
            obj.FaceID = zeros(obj.TotalFaces, 1);
            obj.FaceNormal = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceCenter = zeros(obj.TotalFaces, obj.Dimension);
            obj.FaceArea = zeros(obj.TotalFaces, 1);
            % Vertex Arrays
            obj.VertexCells = cell(obj.TotalVertices, 1);
            obj.VertexFaces = cell(obj.TotalVertices, 1);
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function allocate_more_memory(obj, nverts, ncells, nfaces)
            % Cell Arrays
            if ncells > 0
                new_cells = cell(ncells,1);
                obj.MatID = [obj.MatID;zeros(ncells,1)];
                obj.CellVerts = [obj.CellVerts;cell(ncells,1)];
                obj.CellFaceVerts = [obj.CellFaceVerts;new_cells];
                obj.CellNeighbors = [obj.CellNeighbors;new_cells];
                obj.CellNeighborFaces = [obj.CellNeighborFaces;new_cells];
                obj.CellCenter = [obj.CellCenter;zeros(ncells,obj.Dimension)];
                obj.CellVolume = [obj.CellVolume;zeros(ncells,1)];
                obj.CellSurfaceArea = [obj.CellSurfaceArea;zeros(ncells,1)];
                obj.CellFaces = [obj.CellFaces;new_cells];
                obj.PreviousCell = [obj.PreviousCell;zeros(ncells,1)];
                obj.CellRefinementFlag = [obj.CellRefinementFlag;false(ncells,1)];
                obj.CellRefinementLevel = [obj.CellRefinementLevel;zeros(ncells,1)];
                obj.CellRefinementTop = [obj.CellRefinementTop;zeros(ncells,1)];
                obj.CellRefinementTreeHierarchy = [obj.CellRefinementTreeHierarchy;new_cells];
                obj.RefCellFaces = [obj.RefCellFaces;new_cells];
                obj.RefCellFaceCells = [obj.RefCellFaceCells;new_cells];
                obj.RefCellCornerVerts = [obj.RefCellCornerVerts;new_cells];
                obj.RefCellMidFaceVerts = [obj.RefCellMidFaceVerts;new_cells];
                obj.RefCellHigherLvls = [obj.RefCellHigherLvls;new_cells];
            end
            % Face Arrays
            if nfaces > 0
                new_faces = cell(nfaces,1);
                obj.FaceVerts = [obj.FaceVerts;new_faces];
                obj.FaceCells = [obj.FaceCells;zeros(nfaces,2)];
                obj.OrthogonalProjection = [obj.OrthogonalProjection;zeros(nfaces,2)];
                obj.FaceID = [obj.FaceID;zeros(nfaces,1)];
                obj.FaceNormal = [obj.FaceNormal;zeros(nfaces,obj.Dimension)];
                obj.FaceCenter = [obj.FaceCenter;zeros(nfaces,obj.Dimension)];
                obj.FaceArea = [obj.FaceArea;zeros(nfaces,1)];
            end
            % Vertex Arrays
            if nverts > 0
                new_verts = cell(nverts,1);
                obj.Vertices    = [obj.Vertices;zeros(nverts,obj.Dimension)];
                obj.VertexCells = [obj.VertexCells;new_verts];
                obj.VertexFaces = [obj.VertexFaces;new_verts];
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                              AMR Methods
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function refine_mesh(mesh)
            % ------------------------------------------------------------------------------
            % Get preliminary information
            % ------------------------------------------------------------------------------
            num_old_cells = mesh.TotalCells;
            num_new_cells = sum(mesh.CellRefinementFlag);
            num_new_faces = num_new_cells;
            num_new_verts = num_new_cells;
            mesh.allocate_more_memory(num_new_verts, num_new_cells, num_new_faces);
            % ------------------------------------------------------------------------------
            % Loop through cells and refine cell-by-cell
            % ------------------------------------------------------------------------------
            cc = 0;
            for c=1:num_old_cells
                if mesh.CellRefinementFlag(c)
                    cc = cc + 1;
                    refine_individual_cell(mesh, c, num_old_cells + cc);
                end
            end
            % ------------------------------------------------------------------------------
            % Update Final geometry information
            % ------------------------------------------------------------------------------
            mesh.update_geometry_info_after_modifications();
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function refine_individual_cell( mesh, c, num )
            fvnum = num + 1;
            m = mesh.MatID(c);
            cv1  = mesh.CellVerts{c}(1);  cv2  = mesh.CellVerts{c}(2);
            v1   = mesh.Vertices(cv1,:);  v2   = mesh.Vertices(cv2,:);
            f1   = mesh.CellFaces{c}(1);  f2   = mesh.CellFaces{c}(2);
            fc   = (v1+v2)/2;
            % Face Cells
%             f1c1 = mesh.FaceCells(f1,1); f1c2 = mesh.FaceCells(f1,2);
            f2c1 = mesh.FaceCells(f2,1); f2c2 = mesh.FaceCells(f2,2);
            % ------------------------------------------------------------------------------
            % Make Modifications
            % ------------------------------------------------------------------------------
            mesh.Vertices(fvnum) = fc;
            mesh.MatID(num) = m;
            mesh.FaceVerts{fvnum} = fvnum;
            mesh.FaceID(fvnum) = 0;
            mesh.FaceCells(fvnum,:) = [c,num];
            mesh.CellVerts{c} = [cv1,fvnum];   cv1 = mesh.Vertices(mesh.CellVerts{c});
            mesh.CellVerts{num} = [fvnum,cv2]; cv2 = mesh.Vertices(mesh.CellVerts{num});
            mesh.CellCenter(c,:) = mean(mesh.Vertices(mesh.CellVerts{c}));
            mesh.CellCenter(num,:) = mean(mesh.Vertices(mesh.CellVerts{num}));
            if cv1(2) < cv1(1), fliplr(mesh.CellVerts{c}); end
            if cv2(2) < cv2(1), fliplr(mesh.CellVerts{num}); end
            mesh.CellFaces{c} = [f1,fvnum];
            mesh.CellFaces{num} = [fvnum,f2];
            if f2c1 == c
                mesh.FaceCells(f2,1) = num;
            elseif f2c2 == c
                mesh.FaceCells(f2,2) = num;
            end
            mesh.PreviousCell(c) = c;
            mesh.PreviousCell(num) = c;
            mesh.CellRefinementLevel(c) = mesh.CellRefinementLevel(c) + 1;
            mesh.CellRefinementLevel(num) = mesh.CellRefinementLevel(c);
            mesh.CellRefinementTop(num) = mesh.CellRefinementTop(c);
            % ------------------------------------------------------------------------------
            % Update Refinement Tree - this is actually mostly general for all mesh types
            % ------------------------------------------------------------------------------
            ctop = mesh.CellRefinementTop(c);
            thier = mesh.CellRefinementTreeHierarchy{c}; nthier = length(thier); chier = cell(nthier+1, 1);
            ttree = mesh.CellRefinementTree; tncells = {c,num}; tt = ttree; chier{1} = ttree;
            % Build new hierarchy tree
            for i=1:nthier-1
                ii = i + 1;
                chier{ii,1} = tt{thier(i)};
                tt = tt{thier(i)};
            end
            chier{end,1} = tncells; tt = tncells;
            for i=nthier:-1:1
                ii = thier(i);
                tnew = chier{i};
                tnew{ii} = tt;
                tt = tnew;
            end
            mesh.CellRefinementTree = tt;
            mesh.CellRefinementTreeHierarchy{c} = [thier,1];
            mesh.CellRefinementTreeHierarchy{num} = [thier,2];
            mesh.CellRefinementTreeHierarchyLevel(ctop) = length(mesh.CellRefinementTreeHierarchy{c}) - 1;
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function update_geometry_info_after_modifications(obj)
            % Update All Counts
            obj.TotalVertices = size(obj.Vertices, 1);
            obj.TotalCells = length(obj.CellVerts);
            obj.TotalFaces = length(obj.FaceVerts);
            obj.BoundaryFaces = []; obj.TotalBoundaryFaces = 0;
            obj.InteriorFaces = []; obj.TotalInteriorFaces = 0;
            % Loop through all faces
            for f=1:obj.TotalFaces
                fv = obj.FaceVerts{f};
                fid = obj.FaceID(f);
                if fid == 0
                    obj.TotalInteriorFaces = obj.TotalInteriorFaces + 1;
                    %obj.InteriorFaces = [obj.InteriorFaces;f];
                else
                    obj.TotalBoundaryFaces = obj.TotalBoundaryFaces + 1;
                    %obj.BoundaryFaces = [obj.BoundaryFaces;f];
                end
                fverts = obj.Vertices(fv,:);
                obj.FaceCenter(f,:) = mean(fverts);
                obj.FaceArea(f) = 1;
                fc = obj.FaceCenter(f,:);
                cc = obj.CellCenter(obj.FaceCells(f,1),:);
                obj.FaceNormal(f,:) = (fc-cc)/abs(fc-cc);
            end
            obj.InteriorFaces = zeros(obj.TotalInteriorFaces, 1);
            obj.BoundaryFaces = zeros(obj.TotalBoundaryFaces, 1);
            % Loop through all faces again
            bfcount = 1; ifcount = 1;
            for f=1:obj.TotalFaces
              fid = obj.FaceID(f);
              if fid == 0
                obj.InteriorFaces(ifcount) = f;
                ifcount = ifcount + 1;
              else
                obj.BoundaryFaces(bfcount) = f;
                bfcount = bfcount + 1;
              end
            end
            %obj.TotalBoundaryFaces = length(obj.BoundaryFaces);
            %obj.TotalInteriorFaces = length(obj.InteriorFaces);
            % Loop through all cells
            obj.CellVertexNumbers = [];
            for c=1:obj.TotalCells
                cv = obj.CellVerts{c};
                cverts = obj.Vertices(cv,:);
                cfaces = obj.CellFaces{c};
                obj.CellCenter(c,:) = mean(cverts);
                obj.CellSurfaceArea(c) = sum(obj.FaceArea(cfaces));
                obj.CellVolume(c) = abs(cverts(2) - cverts(1));
            end
            obj.OrthogonalProjection = zeros(obj.TotalFaces,2);
            for f=1:obj.TotalFaces
                fcells = obj.FaceCells(f,:);
                for i=1:2
                    if obj.FaceCells(f,i) == 0, continue; end
                    obj.OrthogonalProjection(f,i) = obj.CellVolume(fcells(i));
                end
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
end