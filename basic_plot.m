
%% test
clear CData CData_vec
% face_list1=getHexorTetFaces(geom1.faces);
% face_list2=getHexorTetFaces(geom2.faces);
CData=geom1.vertices-geom1.vertices_orig;
for counti=1:size(CData,1)
      CData_vec(counti)=norm(CData(counti,:));
end
CData_vec(counti)=norm(CData(counti,:));
% [~,CData_vec]=sort(CData_vec);
% color_set=jet(size(CData,1));
% color_set=color_set(CData_vec,:);
figure()
p1=patch('Faces',geom1.faces,'Vertices',geom1.vertices,'FaceColor','interp','EdgeAlpha',.25,'FaceVertexCData',CData_vec');
hold on
geomtest1=geom1;
geomtest2=geom2;





% for counti=1:size(face_list1,1)
%         current_face=face_list1(counti,:);
%         patch(geom1.vertices(current_face,1),geom1.vertices(current_face,2),geom1.vertices(current_face,3),CData_vec(current_face));
%         hold on
% end
clear CData CData_vec
CData=geom2.vertices-geom2.vertices_orig;
for counti=1:size(CData,1)
      CData_vec(counti)=norm(CData(counti,:));
end
CData_vec(counti)=norm(CData(counti,:));
% [~,CData_vec]=sort(CData_vec);
% color_set=jet(size(CData,1));
% color_set=color_set(CData_vec,:);
p2=patch('Faces',geom2.faces,'Vertices',geom2.vertices,'FaceColor','interp','EdgeAlpha',.3,'FaceVertexCData',CData_vec');
% for counti=1:size(face_list2,1)
%       current_face=face_list2(counti,:);
%       patch(geom2.vertices(current_face,1),geom2.vertices(current_face,2),geom2.vertices(current_face,3),CData_vec(current_face));
%       hold on
% end




colormap jet

% cmap=[0 1 0; 1 .5 0;1 0 0];
% cmap_tot=[];
% c_lin=0:.05:1;
% cmap_tot(:,1)=interp1([0,.5,1]',cmap(:,1),c_lin);
% cmap_tot(:,2)=interp1([0,.5,1]',cmap(:,2),c_lin);
% cmap_tot(:,3)=interp1([0,.5,1]',cmap(:,3),c_lin);
% colormap(cmap_tot)

colorbar
axis equal

hold on
quiver3(geomtest1.vertices_orig(:,1),geomtest1.vertices_orig(:,2),geomtest1.vertices_orig(:,3),...
    (geomtest1.vertices(:,1)-geomtest1.vertices_orig(:,1)),...
    (geomtest1.vertices(:,2)-geomtest1.vertices_orig(:,2)),...
    (geomtest1.vertices(:,3)-geomtest1.vertices_orig(:,3)),'k','AutoScaleFactor',10);

quiver3(geomtest2.vertices_orig(:,1),geomtest2.vertices_orig(:,2),geomtest2.vertices_orig(:,3),...
    (geomtest2.vertices(:,1)-geomtest2.vertices_orig(:,1)),...
    (geomtest2.vertices(:,2)-geomtest2.vertices_orig(:,2)),...
    (geomtest2.vertices(:,3)-geomtest2.vertices_orig(:,3)),'k','AutoScaleFactor',10);


p3=patch('Faces',geomtest1.faces,'Vertices',geomtest1.vertices_orig,'FaceAlpha',.2,'EdgeAlpha',.25);
p4=patch('Faces',geomtest2.faces,'Vertices',geomtest2.vertices_orig,'FaceAlpha',.2,'EdgeAlpha',.25);