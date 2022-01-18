function [faces,vertices]=stlRead2(filename)
        temp=stlread(filename);
        faces=temp.ConnectivityList;
        vertices=temp.Points;
end