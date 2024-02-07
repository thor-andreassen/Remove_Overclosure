function [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure,history_params]=...
    removeOverclosureNODAL(geom1,geom2,params)
    %% Define Gap Threshold

    

    % check values and set defaults
    params=setDefaultParamValue(params,'use_parallel_loops',0);
    params=setDefaultParamValue(params,'relative_gap_weight',0.5);
    params=setDefaultParamValue(params,'smooth_2D_surface',0);
    params=setDefaultParamValue(params,'smoothing',10);
    params=setDefaultParamValue(params,'smoothing_reduction',0.995);
    params=setDefaultParamValue(params,'plot_surf',1);
    params=setDefaultParamValue(params,'stop_tolerance',1E-5);
    params=setDefaultParamValue(params,'geom1_mesh_reduction_factor',0.01);
    params=setDefaultParamValue(params,'geom2_mesh_reduction_factor',0.01);
    params=setDefaultParamValue(params,'scale_reduction_factor',1.005);
    params=setDefaultParamValue(params,'weight_factor',10.0);
    params=setDefaultParamValue(params,'check_original',1);

    % gap to achieve in final meshes in unit of mesh
    desired_gap=params.desired_gap;
    relative_gap_weight=params.relative_gap_weight;
    element_3d_type=params.element_3d_type;
    use_parallel_loops=params.use_parallel_loops;
    plot_surf=params.plot_surf;
    stop_tolerance=abs(params.stop_tolerance);
    check_original=params.check_original;


        % reduction factor initial
        geom1_mesh_reduction_factor=params.geom1_mesh_reduction_factor;
        geom2_mesh_reduction_factor=params.geom2_mesh_reduction_factor;
        scale_percent_factor=params.scale_reduction_factor;


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
    conv_tol=.00001;
    counter=1;
    geom1_error_total=[];
    geom2_error_total=[];
    end_flag=0;


    full_geom1_num_over=[];
    full_geom2_num_over=[];
    full_geom1_max_over=[];
    full_geom2_max_over=[];
    over_total=[];
    if check_original==1
        [full_mesh_params]=calculateFullOverclosure(geom1,geom2,use_parallel_loops,element_3d_type,desired_gap);
        full_geom1_num_over=[full_geom1_num_over,full_mesh_params.geom1_num_over];
        full_geom2_num_over=[full_geom2_num_over,full_mesh_params.geom2_num_over];
        full_geom1_max_over=[full_geom1_max_over,full_mesh_params.geom1_max_over];
        full_geom2_max_over=[full_geom2_max_over,full_mesh_params.geom2_max_over];
    end



    end_flag=0;
    overall_tic=tic();
    overall_time=[];


    while total_error < -stop_tolerance && counter<max_iters && end_flag==0

        %% Reduced Mesh
        rand_ratio=1;
        if element_3d_type(1)
                % element is 3D
                [face_outer_surf,~,~,~,~,geom1_inner_nodes]=get3DElementOuterSurface(geom1.faces,geom1.vertices);
                [geom1.faces_reduce,geom1.vertices_reduce]=renumberFacesAndVertices(face_outer_surf,geom1.vertices);
                temp_rand=randperm(size(geom1_inner_nodes,1));
                temp_rand_val=temp_rand(1:ceil(size(geom1_inner_nodes,1)*rand_ratio));
%                 geom1.vertices_rand=[geom1.vertices_reduce;geom1.vertices(temp_rand_val,:)];
                geom1.vertices_rand=[geom1.vertices];
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
%                 geom2.vertices_rand=[geom2.vertices_reduce;geom2.vertices(temp_rand_val,:)];
                geom2.vertices_rand=geom2.vertices;
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
        scale_val=1;
        geom1_surf_to_geom2_point_distance=geom1_surf_to_geom2_point_distance-desired_gap;
        geom1_surf_to_geom2_vector=geom1_surf_to_geom2_point_points;
        for count_vertex=1:length(geom1_surf_to_geom2_point_distance)
            if geom1_surf_to_geom2_point_distance(count_vertex)>0
                geom1_surf_to_geom2_point_distance(count_vertex)=0;
                multiplier=0;
            else
                multiplier=1;
            end
            geom1_surf_to_geom2_vector(count_vertex,:)=(geom1_surf_to_geom2_point_points(count_vertex,:)...
                -geom2.vertices_rand(count_vertex,:))*multiplier;
            if norm(geom1_surf_to_geom2_vector(count_vertex,:))~=0
                geom1_surf_to_geom2_vector(count_vertex,:)=(geom1_surf_to_geom2_vector(count_vertex,:)/...
                    norm(geom1_surf_to_geom2_vector(count_vertex,:)))*norm(geom1_surf_to_geom2_point_distance(count_vertex))*scale_val;
            end
        end

        geom2_surf_to_geom1_point_distance=geom2_surf_to_geom1_point_distance-desired_gap;
        geom2_surf_to_geom1_vector=geom2_surf_to_geom1_point_points;
        for count_vertex=1:length(geom2_surf_to_geom1_point_distance)
            if geom2_surf_to_geom1_point_distance(count_vertex)>0
                geom2_surf_to_geom1_point_distance(count_vertex)=0;
                multiplier=0;
            else
                multiplier=1;
            end
            geom2_surf_to_geom1_vector(count_vertex,:)=(geom2_surf_to_geom1_point_points(count_vertex,:)...
                -geom1.vertices_rand(count_vertex,:))*multiplier;
            if norm(geom2_surf_to_geom1_vector(count_vertex,:))~=0
                geom2_surf_to_geom1_vector(count_vertex,:)=(geom2_surf_to_geom1_vector(count_vertex,:)/...
                    norm(geom2_surf_to_geom1_vector(count_vertex,:)))*norm(geom2_surf_to_geom1_point_distance(count_vertex))*scale_val;
            end
        end
        
        
        if counter==1
            original_max_overclosure_1=min(geom2_surf_to_geom1_point_distance);
            original_max_overclosure_2=min(geom1_surf_to_geom2_point_distance);
            original_max_overclosure=min([original_max_overclosure_1,original_max_overclosure_2]);
            if original_max_overclosure<-15
                error('gap exceeds 10 mm calculation of over-clsoure');
                disp('the over-closure measure is too large for the current model') 
            end
        end
        
        
        % deform vector is vector needed to deform nodes located from the
        % Surface to the ndoes wheres teh vector is defined. 
%         geom_master_deform_vector=[geom1_surf_to_geom2_vector;-geom2_surf_to_geom1_vector];
%         geom_master_positions=[geom2.vertices_rand;geom1.vertices_rand];


        %% Create arrow deformation plot
        if plot_surf==1
            if counter==1
                def_figure=figure();
            else
                figure(def_figure);
                clf(def_figure);
            end
            subplot(1,2,1)
            plot3(geom2.vertices_rand(:,1),geom2.vertices_rand(:,2),...
                geom2.vertices_rand(:,3),'r.')
            hold on
            if norm(geom1_surf_to_geom2_vector)>0
                quiver3(geom2.vertices_rand(:,1),geom2.vertices_rand(:,2),geom2.vertices_rand(:,3),...
                    geom1_surf_to_geom2_vector(:,1).*geom1_surf_to_geom2_point_distance,...
                    geom1_surf_to_geom2_vector(:,2).*geom1_surf_to_geom2_point_distance,...
                    geom1_surf_to_geom2_vector(:,3).*geom1_surf_to_geom2_point_distance,'k');
            end
            patch('Faces',geom1.faces_reduce,'Vertices',geom1.vertices_reduce,'FaceAlpha',.3,'EdgeAlpha',.3,'FaceColor','b');
            axis equal

            subplot(1,2,2)
            plot3(geom1.vertices_rand(:,1),geom1.vertices_rand(:,2),...
                geom1.vertices_rand(:,3),'r.')
            hold on
            if norm(geom2_surf_to_geom1_vector)>0
                quiver3(geom1.vertices_rand(:,1),geom1.vertices_rand(:,2),geom1.vertices_rand(:,3),...
                    geom2_surf_to_geom1_vector(:,1).*geom2_surf_to_geom1_point_distance,...
                    geom2_surf_to_geom1_vector(:,2).*geom2_surf_to_geom1_point_distance,...
                    geom2_surf_to_geom1_vector(:,3).*geom2_surf_to_geom1_point_distance,'k');
            end
            patch('Faces',geom2.faces_reduce,'Vertices',geom2.vertices_reduce,'FaceAlpha',.3,'EdgeAlpha',.3,'FaceColor','b');
            axis equal
            pause(.001);
        end


        %% Apply deformations
        geom2.vertices=geom2.vertices+geom1_surf_to_geom2_vector;
        geom1.vertices=geom1.vertices+geom2_surf_to_geom1_vector;
% % %         %% Display Original Min Gap
% % %         geom2_error=min(geom1_surf_to_geom2_point_distance);
% % %         geom1_error=min(geom2_surf_to_geom1_point_distance);
% % %         geom1_error_total=[geom1_error_total,geom1_error];
% % %         geom2_error_total=[geom2_error_total,geom2_error];
% % %         
% % %         
% % %         if counter>20
% % %             if geom1_error_total(end)~=0
% % %                 conv_1=abs((geom1_error_total(end)-geom1_error_total(end-1))/geom1_error_total(end));
% % %             else
% % %                 conv_1=0;
% % %             end
% % %             if geom2_error_total(end)~=0
% % %                 conv_2=abs((geom2_error_total(end)-geom2_error_total(end-1))/geom2_error_total(end));
% % %             else
% % %                 conv_2=0;
% % %             end
% % %             
% % %             
% % %            max_conv=max([conv_1,conv_2]);
% % %            if max_conv<=conv_tol
% % %                end_flag=1;
% % %            else
% % %                end_flag=0;
% % %            end
% % %         end
% % %         table(geom2_error,geom1_error)
% % %         total_error=min([geom2_error,geom1_error])

%% Display Original Min Gap
        geom2_error=min(geom1_surf_to_geom2_point_distance);
        geom1_error=min(geom2_surf_to_geom1_point_distance);
        geom1_error_mean=mean(geom2_surf_to_geom1_point_distance);
        geom2_error_mean=mean(geom1_surf_to_geom2_point_distance);
        geom1_error_total=[geom1_error_total,geom1_error];
        geom2_error_total=[geom2_error_total,geom2_error];
        geom_error_mean=log10(abs(min([geom1_error_mean,geom2_error_mean])));

        max_deform=log10(norm([geom1_surf_to_geom2_vector;geom2_surf_to_geom1_vector]));
        max_adjustment_vec=vecnorm([geom1_surf_to_geom2_vector;geom2_surf_to_geom1_vector],2,2);
        current_num_overclosures = log10(length(nonzeros(max_adjustment_vec(max_adjustment_vec>0)))/length(max_adjustment_vec));
        current_over=length(nonzeros(max_adjustment_vec(max_adjustment_vec>0)));
        over_total=[over_total,current_over];

        history_params.geom1_error_total=geom1_error_total;
        history_params.geom2_error_total=geom2_error_total;
        history_params.over_total=over_total;
        if counter>150
            if geom1_error_total(end)~=0
                conv_1=abs((geom1_error_total(end)-geom1_error_total(end-1))/geom1_error_total(end));
            else
                conv_1=0;
            end
            if geom2_error_total(end)~=0
                conv_2=abs((geom2_error_total(end)-geom2_error_total(end-1))/geom2_error_total(end));
            else
                conv_2=0;
            end


            max_conv=max([conv_1,conv_2]);
            if max_conv<=conv_tol
                end_flag=1;
            else
                end_flag=0;
            end
        end
        table(geom2_error,geom1_error,geom1_error_mean,geom2_error_mean,current_over)
        total_error=min([geom2_error,geom1_error])
        if total_error>=-.001
            end_flag=1;
        end


        if counter==1
            conv_figure=figure();
        else
            figure(conv_figure);
            clf(conv_figure);
        end
        subplot(2,2,1);
        total_error_hist(counter)=log10(abs(total_error));
        plot(total_error_hist);
        xlabel('Iteration');
        ylabel('Log MAX Overclosure');
        title('Maximum Overclosure');

        subplot(2,2,2);
        mean_conv_hist(counter)=geom_error_mean;
        plot(mean_conv_hist);
        xlabel('Iteration');
        ylabel('Log Mean Overclosure');
        title('Mean Overclosure');

        subplot(2,2,3);
        max_deform_hist(counter)=max_deform;
        plot(max_deform_hist);
        xlabel('Iteration');
        ylabel('Log Max Adjustment');
        title('Max Adjustment');

        subplot(2,2,4);
        num_over(counter)=current_num_overclosures;
        plot(num_over);
        xlabel('Iteration');
        ylabel('log Fraction of Overclosures');
        title('Fraction of Overclosures');



        scale_val=scale_val*.99;
        counter=counter+1;

    if check_original==1
            [full_mesh_params]=calculateFullOverclosure(geom1,geom2,use_parallel_loops,element_3d_type,desired_gap);
            full_geom1_num_over=[full_geom1_num_over,full_mesh_params.geom1_num_over];
            full_geom2_num_over=[full_geom2_num_over,full_mesh_params.geom2_num_over];
            full_geom1_max_over=[full_geom1_max_over,full_mesh_params.geom1_max_over];
            full_geom2_max_over=[full_geom2_max_over,full_mesh_params.geom2_max_over];
            history_params.full_geom1_num_over=full_geom1_num_over;
            history_params.full_geom2_num_over=full_geom2_num_over;
            history_params.full_geom1_max_over=full_geom1_max_over;
            history_params.full_geom2_max_over=full_geom2_max_over;
        end
        history_params.overall_time=overall_time;


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