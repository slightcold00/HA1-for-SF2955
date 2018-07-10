%calculate euclidean
function myOutput = cal(l,x,pos_vec)
    part_final =  [x(1,:); x(4,:)]' -pos_vec(:,l)';
    len = length(x(1,:));
    eu_dist = zeros(len,1);
    for i=1:len
        eu_dist(i,1) = sqrt((part_final(i,1)^2)+((part_final(i,2))^2)); 
        myOutput = transpose(eu_dist);
    end
end