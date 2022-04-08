%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate distance for horizontal eddy calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear distu_2D distv_2D
for i=1:NY-1
for j=1:NX-1
        distu_2D(i,j) = (lon(i,j) - lon(i,j+1)).^2 + (lat(i,j) - lat(i,j+1)).^2;
        distv_2D(i,j) = (lon(i,j) - lon(i+1,j)).^2 + (lat(i,j) - lat(i+1,j)).^2;
end
end

for i=1:NY-1
for j=1:NX-1
        distu_2D(i,j) = sw_dist([lat(i,j) lat(i,j+1)] , [lon(i,j) lon(i,j+1)] ,'km').*1e3 ;
        distv_2D(i,j) = sw_dist([lat(i,j) lat(i+1,j)] , [lon(i,j) lon(i+1,j)] , 'km').*1e3 ;
end
end

du1 = nan(NY,NX);
dv1 = nan(NY,NX);
du1(2:end,2:end) = distu_2D ;
dv1(2:end,2:end) = distv_2D ;

du2 = repmat(du1,1,1,NZ);
du2 = permute(du2, [3 1 2]);
dv2 = repmat(dv1,1,1,NZ);
dv2 = permute(dv2, [3 1 2]);

