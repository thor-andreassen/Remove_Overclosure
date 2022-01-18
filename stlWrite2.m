function stlWrite2(filename,faces,vertices)
        temp=triangulation(faces,vertices);
        stlwrite(temp,filename);
end