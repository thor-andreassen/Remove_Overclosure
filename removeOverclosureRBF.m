function [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure]=...
    removeOverclosureRBF(geom1,geom2,params)
    %% Define Gap Threshold

        % gap to achieve in final meshes in unit of mesh
        desired_gap=params.desired_gap;
        %1 = fixed geom 1 moving geom 2.
        % 0 = moving geom 1, fixed geom 2.
        relative_gap_weight=params.relative_gap_weight;
        element_3d_type=params.element_3d_type;
        use_parallel_loops=params.use_parallel_loops;


        smoothing_improve=params.smoothing_improve;
        plot_surf=params.plot_surf;
        % good parameters for 2D and 3D
    %     smoothing=50;
    %     rbf_iterations=300;
        smoothing=params.smoothing;
        rbf_iterations=params.rbf_iterations;


        % reduction factor initial
        geom1_mesh_reduction_factor=params.geom1_mesh_reduction_factor;
        geom2_mesh_reduction_factor=params.geom2_mesh_reduction_factor;
        scale_percent_factor=params.scale_percent_factor;

        % good parameters for 2D and 2D
    %         smoothing=1000;
    %         rbf_iterations=1000;
    %         smoothing=1000;
    %         rbf_iterations=400;
    %% load initial geometries
    %     load stls
    %     [geom2.faces, geom2.vertices]=stlRead2('VHFL_Left_Bone_Femur.stl');
    %     [geom1.faces, geom1.vertices]=stlRead2('VHFL_Muscle_VastusIntermedius.stl');
    % 
    %     save('muscle_geom_orig.mat');      

    %% main
    % input to the function is two variables called geom1 and geom 2.
    %These are structures that have the folloiwng variables:
            % geom.faces: A connectivity list representing the set
                    % of elements. Each row represents a different element,
                    % with the row number corresponding to the assumed row.
                    % The elements can be of the following types and
                    % sizes:
                            % tri: (E x 3) triangular mesh with
                                    % 3 nodes per face of a triangular
                                    % mesh (3 nodes, 3 edges, 1
                                    % face)
                            % quad: (E x 4) quadrilateral mesh
                                    % with 4 nodes per face of a
                                    % quadrilateral mesh.
                                    % (4 nodes, 4 edges, 1
                                    % face)
                            % tetrahedral: (E x 4) tetrahedral
                                    % mesh with 4 nodes per element.
                                    % (4 nodes, 6 edges, 4
                                    % faces)
                            % hexahedral: (E x 8) hexahedral
                                    % mesh with 8 nodes per element.
                                    % (8 nodes, 12 edges, 6
                                    % faces)
            % geom.vertices: (m x 3) A list of nodes as vertices of the
                    % elements with the row corresponding to the assumed node number
                    % corresponding to the list of nodes in the elements. The values
                    % represent the location in cartesian x, y, z coordinate space.

    % the other input is a binary array of [1 x 2] with a definitions of
            % whether the elements inputted to each geometry are 3d or 2d respectively.

    total_error=-Inf;
    max_iters=150;
    counter=1;
    while total_error < -.1 && counter<max_iters

        %% Reduced Mesh
        rand_ratio=.75;
        if element_3d_type(1)
                % element is 3D
                [face_outer_surf,~,~,~,~,geom1_inner_nodes]=get3DElementOuterSurface(geom1.faces,geom1.vertices);
                [geom1.faces_reduce,geom1.vertices_reduce]=renumberFacesAndVertices(face_outer_surf,geom1.vertices);
                temp_rand=randperm(size(geom1_inner_nodes,1));
                temp_rand_val=temp_rand(1:ceil(size(geom1_inner_nodes,1)*rand_ratio));
                geom1.vertices_rand=[geom1.vertices_reduce;geom1.vertices(temp_rand_val,:)];
                if size(geom1.faces_reduce,2)==4
                        geom1_reduce_type_Q4=1;
                else
                        geom1_reduce_type_Q4=0;
                end
        else
                % element is 2D
                if size(geom1.faces,2)==3
                    
                        geom1_reduce_type_Q4=0;
                        % element is a tri
    %                     geom1_mesh_reduction_factor=750;
    %                     geom1_mesh_reduction_factor=2000;
                        if counter==1
                            try
                                if smoothing_improve==1 && geom1_reduce_type_Q4==0 && element_3d_type(1)==0 && relative_gap_weight~=1
                                    geom1.vertices=improveTriMeshQuality(geom1.faces,geom1.vertices,2,2,.001);
                                end
                            catch
                                disp('geom1 mesh improvement failed');
                            end
                        end
                       geom1_mesh_reduction_factor=scaleInputReductionFactor(geom1_mesh_reduction_factor,scale_percent_factor);
                       if  geom1_mesh_reduction_factor<1
                           temp=reducepatch(geom1.faces,geom1.vertices,geom1_mesh_reduction_factor);
                            geom1.faces_reduce=temp.faces;
                            geom1.vertices_reduce=temp.vertices;
                       else
                           geom1.faces_reduce=geom1.faces;
                            geom1.vertices_reduce=geom1.vertices;
                       end
                else
                        % element is a quad
                        geom1_reduce_type_Q4=1;
                end
                geom1.vertices_rand=geom1.vertices_reduce;
        end

       if element_3d_type(2)
               % element is 3D
               [face_outer_surf,~,~,~,~,geom2_inner_nodes]=get3DElementOuterSurface(geom2.faces,geom2.vertices);
                [geom2.faces_reduce,geom2.vertices_reduce]=renumberFacesAndVertices(face_outer_surf,geom2.vertices);
                temp_rand=randperm(size(geom2_inner_nodes,1));
                temp_rand_val=temp_rand(1:ceil(size(geom2_inner_nodes,1)*rand_ratio));
                geom2.vertices_rand=[geom2.vertices_reduce;geom2.vertices(temp_rand_val,:)];
                if size(geom2.faces_reduce,2)==4
                        geom2_reduce_type_Q4=1;
                else
                        geom2_reduce_type_Q4=0;
                end
       else
               % element is 2D
               if size(geom2.faces,2)==3
                   
                       geom2_reduce_type_Q4=0;
                        % element is a tri
    %                     geom2_mesh_reduction_factor=750;
    %                    geom2_mesh_reduction_factor=2000;
                        if counter==1
                            try
                                if smoothing_improve==1 && geom2_reduce_type_Q4==0 && element_3d_type(2)==0 && relative_gap_weight~=0
                                    geom2.vertices=improveTriMeshQuality(geom2.faces,geom2.vertices,2,2,.001);
                                end
                            catch
                                disp('geom2 mesh improvement failed');
                            end
                        end
                       geom2_mesh_reduction_factor=scaleInputReductionFactor(geom2_mesh_reduction_factor,scale_percent_factor);
                       if geom2_mesh_reduction_factor<1
                           temp=reducepatch(geom2.faces,geom2.vertices,geom2_mesh_reduction_factor);
                           geom2.faces_reduce=temp.faces;
                           geom2.vertices_reduce=temp.vertices;
                       else
                           geom2.faces_reduce=geom2.faces;
                           geom2.vertices_reduce=geom2.vertices;
                       end
               else
                       % element is a quad
                       geom2_reduce_type_Q4=1;
               end
               geom2.vertices_rand=geom2.vertices_reduce;
       end


        %% initial plots
    %     geom2_fig=figure()
    %     geom2_patch_orig=patch('Faces',geom2.faces,'Vertices',geom2.vertices,'FaceColor','r','FaceAlpha',.5)
    %     hold on
    %     geom2_patch_reduce=patch('Faces',geom2.faces_reduce,'Vertices',geom2.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    %     
    %     geom1_fig=figure()
    %     geom1_patch_orig=patch('Faces',geom1.faces,'Vertices',geom1.vertices,'FaceColor','r','FaceAlpha',.5)
    %     hold on
    %     geom1_patch_reduce=patch('Faces',geom1.faces_reduce,'Vertices',geom1.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    %     
    %     both_geom_fig=figure()
    %     geom2_patch_orig=patch('Faces',geom2.faces,'Vertices',geom2.vertices,'FaceColor','r','FaceAlpha',.5)
    %     hold on
    %     geom2_patch_reduce=patch('Faces',geom2.faces_reduce,'Vertices',geom2.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    %     geom1_patch_orig=patch('Faces',geom1.faces,'Vertices',geom1.vertices,'FaceColor','r','FaceAlpha',.5)
    %     hold on
    %     geom1_patch_reduce=patch('Faces',geom1.faces_reduce,'Vertices',geom1.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    %     
    
    %% Initial geometry smoothing

    
    
        %% determine initial overclosure distances
        if geom1_reduce_type_Q4==0
            % mesh geometry is tri
            tri_timer=tic();
            [geom1_surf_to_geom2_point_distance, geom1_surf_to_geom2_point_points] = point2trimesh('Faces',geom1.faces_reduce,...
                'Vertices',geom1.vertices_reduce,'QueryPoints',geom2.vertices_rand,'Algorithm','parallel');
%             toc(tri_timer)
        elseif geom1_reduce_type_Q4==1
            % mesh geometry is quad
            [geom1_surf_to_geom2_point_points,geom1_surf_to_geom2_point_distance]=...
                    getPointToQ4MeshApproximate(geom1.faces_reduce,geom1.vertices_reduce,geom2.vertices_rand,use_parallel_loops);
        end

        if geom2_reduce_type_Q4==0
            % mesh geometry is tri
            tri_timer=tic();
            [geom2_surf_to_geom1_point_distance, geom2_surf_to_geom1_point_points] = point2trimesh('Faces',geom2.faces_reduce,...
                'Vertices',geom2.vertices_reduce,'QueryPoints',geom1.vertices_rand,'Algorithm','parallel');

%             toc(tri_timer)
        elseif geom2_reduce_type_Q4==1
           % mesh geometry is quad
           [geom2_surf_to_geom1_point_points,geom2_surf_to_geom1_point_distance]=...
                    getPointToQ4MeshApproximate(geom2.faces_reduce,geom2.vertices_reduce,geom1.vertices_rand,use_parallel_loops);
        end


        %% Determine Gaps

        geom1_surf_to_geom2_point_distance_dir=sign(geom1_surf_to_geom2_point_distance);
        geom1_surf_to_geom2_point_distance=geom1_surf_to_geom2_point_distance-desired_gap;
        geom1_surf_to_geom2_vector=geom1_surf_to_geom2_point_points;
        for count_vertex=1:length(geom1_surf_to_geom2_point_distance)
            if geom1_surf_to_geom2_point_distance(count_vertex)>0
                geom1_surf_to_geom2_point_distance(count_vertex)=0;
                multiplier=0;
            else
                multiplier=1;
            end
            geom1_surf_to_geom2_vector(count_vertex,:)=-geom1_surf_to_geom2_point_distance_dir(count_vertex)*(geom1_surf_to_geom2_point_points(count_vertex,:)...
                -geom2.vertices_rand(count_vertex,:))*multiplier;
        end

        geom2_surf_to_geom1_point_distance_dir=sign(geom2_surf_to_geom1_point_distance);
        geom2_surf_to_geom1_point_distance=geom2_surf_to_geom1_point_distance-desired_gap;
        geom2_surf_to_geom1_vector=geom2_surf_to_geom1_point_points;
        for count_vertex=1:length(geom2_surf_to_geom1_point_distance)
            if geom2_surf_to_geom1_point_distance(count_vertex)>0
                geom2_surf_to_geom1_point_distance(count_vertex)=0;
                multiplier=0;
            else
                multiplier=1;
            end
            geom2_surf_to_geom1_vector(count_vertex,:)=-geom2_surf_to_geom1_point_distance_dir(count_vertex)*(geom2_surf_to_geom1_point_points(count_vertex,:)...
                -geom1.vertices_rand(count_vertex,:))*multiplier;
        end
        
        
        if counter==1
            original_max_overclosure_1=min(geom2_surf_to_geom1_point_distance);
            original_max_overclosure_2=min(geom1_surf_to_geom2_point_distance);
            original_max_overclosure=min([original_max_overclosure_1,original_max_overclosure_2]);
            if original_max_overclosure<-10
                error('gap exceeds 10 mm calculation of over-clsoure');
                disp('the over-closure measure is too large for the current model') 
            end
        end
        
        
        % deform vector is vector needed to deform nodes located from the
        % Surface to the ndoes wheres teh vector is defined. 
        geom_master_deform_vector=[geom1_surf_to_geom2_vector;-geom2_surf_to_geom1_vector];
        geom_master_positions=[geom2.vertices_rand;geom1.vertices_rand];

    %     geom1_surf_to_geom2_point_distance_dir=sign(geom1_surf_to_geom2_point_distance);
    %     geom1_surf_to_geom2_point_distance=geom1_surf_to_geom2_point_distance-desired_gap;
    %     geom1_surf_to_geom2_vector=geom1_surf_to_geom2_point_points;
    %     for count_vertex=1:length(geom1_surf_to_geom2_point_distance)
    %         if geom1_surf_to_geom2_point_distance(count_vertex)>0
    %             geom1_surf_to_geom2_point_distance(count_vertex)=0;
    %             multiplier=0;
    %         else
    %             multiplier=1;
    %         end
    %         geom1_surf_to_geom2_vector(count_vertex,:)=-geom1_surf_to_geom2_point_distance_dir(count_vertex)*(geom1_surf_to_geom2_point_points(count_vertex,:)...
    %             -geom2.vertices_reduce(count_vertex,:))*(1-relative_gap_weight)*multiplier;
    %     end
    % 
    % 
    % 
    %     geom2_surf_to_geom1_point_distance_dir=sign(geom2_surf_to_geom1_point_distance);
    %     geom2_surf_to_geom1_point_distance=geom2_surf_to_geom1_point_distance-desired_gap;
    %     geom2_surf_to_geom1_vector=geom2_surf_to_geom1_point_points;
    %     for count_vertex=1:length(geom2_surf_to_geom1_point_distance)
    %         if geom2_surf_to_geom1_point_distance(count_vertex)>0
    %             geom2_surf_to_geom1_point_distance(count_vertex)=0;
    %             multiplier=0;
    %         else
    %             multiplier=1;
    %         end
    %         geom2_surf_to_geom1_vector(count_vertex,:)=-geom2_surf_to_geom1_point_distance_dir(count_vertex)*(geom2_surf_to_geom1_point_points(count_vertex,:)...
    %             -geom1.vertices_reduce(count_vertex,:))*relative_gap_weight*multiplier;
    %     end

        %% Create arrow deformation plot
        if plot_surf==1
            geom2_deform_reduce_fig=figure();
            plot3(geom2.vertices_rand(:,1),geom2.vertices_rand(:,2),...
                geom2.vertices_rand(:,3),'ro')
            hold on
            if norm(geom1_surf_to_geom2_vector)>0
        %         arrow3(geom2.vertices_rand,geom1_surf_to_geom2_vector+geom2.vertices_rand)
                quiver3(geom2.vertices_rand(:,1),geom2.vertices_rand(:,2),geom2.vertices_rand(:,3),...
                    geom1_surf_to_geom2_vector(:,1).*geom1_surf_to_geom2_point_distance,...
                    geom1_surf_to_geom2_vector(:,2).*geom1_surf_to_geom2_point_distance,...
                    geom1_surf_to_geom2_vector(:,3).*geom1_surf_to_geom2_point_distance,'k');
            end
            patch('Faces',geom1.faces_reduce,'Vertices',geom1.vertices_reduce,'FaceAlpha',.3,'EdgeAlpha',.3,'FaceColor','b');
            axis equal

            geom1_deform_reduce_fig=figure();
            plot3(geom1.vertices_rand(:,1),geom1.vertices_rand(:,2),...
                geom1.vertices_rand(:,3),'ro')
            hold on
            if norm(geom2_surf_to_geom1_vector)>0
        %         arrow3(geom1.vertices_rand,geom2_surf_to_geom1_vector+geom1.vertices_rand)
                quiver3(geom1.vertices_rand(:,1),geom1.vertices_rand(:,2),geom1.vertices_rand(:,3),...
                    -geom2_surf_to_geom1_vector(:,1).*geom2_surf_to_geom1_point_distance,...
                    -geom2_surf_to_geom1_vector(:,2).*geom2_surf_to_geom1_point_distance,...
                    -geom2_surf_to_geom1_vector(:,3).*geom2_surf_to_geom1_point_distance,'k');
            end
            patch('Faces',geom2.faces_reduce,'Vertices',geom2.vertices_reduce,'FaceAlpha',.3,'EdgeAlpha',.3,'FaceColor','b');
            axis equal
        end
        %% create Radial Basis Approximation

%         rbf_timer=tic();



%     geom_deform_vec_rbf=newrb(geom_master_positions',geom_master_deform_vector',...
%             1E-6,smoothing,rbf_iterations);
    smoothing=smoothing*.9;
    geom_deform_vec_rbf=newgrnn(geom_master_positions',geom_master_deform_vector',smoothing);
    
    
    %     geom_deform_vec_rbf=newrbe(geom_master_positions',geom_master_deform_vector',1000);



        %% Determine Original Deformations
        if use_parallel_loops
                geom2_deform_orig_vec=sim(geom_deform_vec_rbf,geom2.vertices','useParallel','yes');
                geom2_deform_orig_vec=geom2_deform_orig_vec'*(relative_gap_weight);
                geom1_deform_orig_vec=sim(geom_deform_vec_rbf,geom1.vertices','useParallel','yes');
                geom1_deform_orig_vec=-geom1_deform_orig_vec'*(1-relative_gap_weight);

        else
                geom2_deform_orig_vec=sim(geom_deform_vec_rbf,geom2.vertices');
                geom2_deform_orig_vec=geom2_deform_orig_vec'*(relative_gap_weight);
                geom1_deform_orig_vec=sim(geom_deform_vec_rbf,geom1.vertices');
                geom1_deform_orig_vec=-geom1_deform_orig_vec'*(1-relative_gap_weight);
        end
        %% Plot original deformations
    %     geom2_deform_orig_fig=figure();
    %     for count_face=1:size(geom2.faces,1)
    %         nodel=geom2.faces(count_face,:);
    %         temp_vec=geom2_deform_orig_vec(nodel,:);
    %         temp_vec_mag=zeros(3,1);
    %         for count_vec=1:3
    %             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
    %         end
    %         patch(geom2.vertices(nodel,1),geom2.vertices(nodel,2),...
    %             geom2.vertices(nodel,3),temp_vec_mag);
    %         hold on
    %     end
    %     colorbar
    %     colormap jet
    %     
    %     
    %     geom1_deform_orig_fig=figure();
    %     for count_face=1:size(geom1.faces,1)
    %         nodel=geom1.faces(count_face,:);
    %         temp_vec=geom1_deform_orig_vec(nodel,:);
    %         temp_vec_mag=zeros(3,1);
    %         for count_vec=1:3
    %             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
    %         end
    %         patch(geom1.vertices(nodel,1),geom1.vertices(nodel,2),...
    %             geom1.vertices(nodel,3),temp_vec_mag);
    %         hold on
    %     end
    %     colorbar
    %     colormap jet


        %% Apply deformations
    %     deforgeom22_fig=figure()
    %     plot3(geom2.vertices(:,1),geom2.vertices(:,2),...
    %         geom2.vertices(:,3),'ro');
    %     hold on
        geom2.vertices=geom2.vertices+geom2_deform_orig_vec;
    %     plot3(geom2.vertices(:,1),geom2.vertices(:,2),...
    %         geom2.vertices(:,3),'bo');
    %     if norm(geom1_surf_to_geom2_vector)>0
    %         arrow3(geom2.vertices_reduce,geom1_surf_to_geom2_vector+geom2.vertices_reduce)
    %     end
    % 
    %     deforgeom21_fig=figure()
    %     plot3(geom1.vertices(:,1),geom1.vertices(:,2),...
    %         geom1.vertices(:,3),'ro');
    %     hold on
        geom1.vertices=geom1.vertices+geom1_deform_orig_vec;
    %     plot3(geom1.vertices(:,1),geom1.vertices(:,2),...
    %         geom1.vertices(:,3),'bo');
    %     if norm(geom2_surf_to_geom1_vector)>0
    %         arrow3(geom1.vertices_reduce,geom2_surf_to_geom1_vector+geom1.vertices_reduce+.0001)
    %     end
        %% Display Original Min Gap
        geom2_error=min(geom1_surf_to_geom2_point_distance);
        geom1_error=min(geom2_surf_to_geom1_point_distance);
        table(geom2_error,geom1_error)
        total_error=min([geom2_error,geom1_error]);
        %% Save new stls
        counter=counter+1;
%         try
%             if mod(counter,10)==0 && smoothing_improve==1 && geom1_reduce_type_Q4==0 && element_3d_type(1)==0 && relative_gap_weight~=1
%                 geom1.vertices=improveTriMeshQuality(geom1.faces,geom1.vertices,2,1,.001);
%             end
%         catch
%             disp('geom1 mesh improvement failed');
%         end
% 
% 
%         try
%             if mod(counter,10)==0 && smoothing_improve==1 && geom2_reduce_type_Q4==0 && element_3d_type(2)==0 && relative_gap_weight~=0
%                 geom2.vertices=improveTriMeshQuality(geom2.faces,geom2.vertices,2,1,.001);
%             end
%         catch
%             disp('geom1 mesh improvement failed');
%         end
    end


%     try
%         if smoothing_improve==1 && geom1_reduce_type_Q4==0 && element_3d_type(1)==0 && relative_gap_weight~=1
%             geom1.vertices=improveTriMeshQuality(geom1.faces,geom1.vertices,2,2,.001);
%         end
%     catch
%         disp('geom1 mesh improvement failed');
%     end
% 
% 
%     try
%         if smoothing_improve==1 && geom2_reduce_type_Q4==0 && element_3d_type(2)==0 && relative_gap_weight~=0
%             geom2.vertices=improveTriMeshQuality(geom2.faces,geom2.vertices,2,2,.001);
%         end
%     catch
%         disp('geom1 mesh improvement failed');
%     end
    geom1_new=geom1;
    geom2_new=geom2;
%     if counter>10
%         error('max iterations exceeded');
%     end
end