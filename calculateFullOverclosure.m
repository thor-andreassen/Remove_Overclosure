function [full_mesh_params]=calculateFullOverclosure(geom1,geom2,use_parallel_loops,element_3d_type,desired_gap)
    counter=0;    
    %% Reduced Mesh
    rand_ratio=1;
    if element_3d_type(1)
        % element is 3D
        [face_outer_surf,~,~,~,~,geom1_inner_nodes]=get3DElementOuterSurface(geom1.elems,geom1.vertices);
        [geom1.elems_reduce,geom1.vertices_reduce]=renumberFacesAndVertices(face_outer_surf,geom1.vertices);
        temp_rand=randperm(size(geom1_inner_nodes,1));
        temp_rand_val=temp_rand(1:ceil(size(geom1_inner_nodes,1)*rand_ratio));
        geom1.vertices_rand=[geom1.vertices_reduce;geom1_inner_nodes(temp_rand_val,:)];
        if size(geom1.elems_reduce,2)==4
            geom1_reduce_type_Q4=1;
        else
            geom1_reduce_type_Q4=0;
        end
    else
        % element is 2D
        if size(geom1.elems,2)==3
    
            geom1_reduce_type_Q4=0;
            if counter==1
                try
                    if smooth_2D_surface==1 && geom1_reduce_type_Q4==0 && element_3d_type(1)==0 && relative_gap_weight~=1
                        geom1.vertices=improveTriMeshQuality(geom1.elems,geom1.vertices,2,2,.001);
                    end
                catch
                    disp('geom1 mesh improvement failed');
                end
            end
            geom1_mesh_reduction_factor=1;
            if  geom1_mesh_reduction_factor<1
                temp=reducepatch(geom1.elems,geom1.vertices,geom1_mesh_reduction_factor);
                geom1.elems_reduce=temp.elems;
                geom1.vertices_reduce=temp.vertices;
            else
                geom1.elems_reduce=geom1.elems;
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
        [face_outer_surf,~,~,~,~,geom2_inner_nodes]=get3DElementOuterSurface(geom2.elems,geom2.vertices);
        [geom2.elems_reduce,geom2.vertices_reduce]=renumberFacesAndVertices(face_outer_surf,geom2.vertices);
        temp_rand=randperm(size(geom2_inner_nodes,1));
        temp_rand_val=temp_rand(1:ceil(size(geom2_inner_nodes,1)*rand_ratio));
        geom2.vertices_rand=[geom2.vertices_reduce;geom2_inner_nodes(temp_rand_val,:)];
        if size(geom2.elems_reduce,2)==4
            geom2_reduce_type_Q4=1;
        else
            geom2_reduce_type_Q4=0;
        end
    else
        % element is 2D
        if size(geom2.elems,2)==3
            geom2_reduce_type_Q4=0;
            if counter==1
                try
                    if smooth_2D_surface==1 && geom2_reduce_type_Q4==0 && element_3d_type(2)==0 && relative_gap_weight~=0
                        geom2.vertices=improveTriMeshQuality(geom2.elems,geom2.vertices,2,2,.001);
                    end
                catch
                    disp('geom2 mesh improvement failed');
                end
            end
            geom2_mesh_reduction_factor=1;
            if geom2_mesh_reduction_factor<1
                temp=reducepatch(geom2.elems,geom2.vertices,geom2_mesh_reduction_factor);
                geom2.elems_reduce=temp.elems;
                geom2.vertices_reduce=temp.vertices;
            else
                geom2.elems_reduce=geom2.elems;
                geom2.vertices_reduce=geom2.vertices;
            end
        else
            % element is a quad
            geom2_reduce_type_Q4=1;
        end
        geom2.vertices_rand=geom2.vertices_reduce;
    end
    
    %% determine initial overclosure distances
    if geom1_reduce_type_Q4==0
        % mesh geometry is tri
        tri_timer=tic();
        [geom1_surf_to_geom2_point_distance, geom1_surf_to_geom2_point_points] = point2trimesh('Faces',geom1.elems_reduce,...
            'Vertices',geom1.vertices_reduce,'QueryPoints',geom2.vertices_rand,'Algorithm','parallel');
        %             toc(tri_timer)
    elseif geom1_reduce_type_Q4==1
        % mesh geometry is quad
        [geom1_surf_to_geom2_point_points,geom1_surf_to_geom2_point_distance]=...
            getPointToQ4MeshApproximate(geom1.elems_reduce,geom1.vertices_reduce,geom2.vertices_rand,use_parallel_loops);
    
    end
    
    
    
    
    if geom2_reduce_type_Q4==0
        % mesh geometry is tri
        tri_timer=tic();
        [geom2_surf_to_geom1_point_distance, geom2_surf_to_geom1_point_points] = point2trimesh('Faces',geom2.elems_reduce,...
            'Vertices',geom2.vertices_reduce,'QueryPoints',geom1.vertices_rand,'Algorithm','parallel');
    
        %             toc(tri_timer)
    elseif geom2_reduce_type_Q4==1
        % mesh geometry is quad
        [geom2_surf_to_geom1_point_points,geom2_surf_to_geom1_point_distance]=...
            getPointToQ4MeshApproximate(geom2.elems_reduce,geom2.vertices_reduce,geom1.vertices_rand,use_parallel_loops);
    end
    
    %% Determine Gaps
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
                norm(geom1_surf_to_geom2_vector(count_vertex,:)))*norm(geom1_surf_to_geom2_point_distance(count_vertex));
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
                norm(geom2_surf_to_geom1_vector(count_vertex,:)))*norm(geom2_surf_to_geom1_point_distance(count_vertex));
        end
    end
    
    full_mesh_params.geom1_num_over=length(nonzeros(geom2_surf_to_geom1_point_distance(geom2_surf_to_geom1_point_distance<0)));
    full_mesh_params.geom2_num_over=length(nonzeros(geom1_surf_to_geom2_point_distance(geom1_surf_to_geom2_point_distance<0)));
    full_mesh_params.geom1_max_over=max(abs(geom2_surf_to_geom1_point_distance));
    full_mesh_params.geom2_max_over=max(abs(geom1_surf_to_geom2_point_distance));

end