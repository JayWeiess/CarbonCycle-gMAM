function [ vec,residue,k] = findp( x1,x2,delta,res )
%Find the equal-distance(delta) points on the segement between vector x1
% and x2, and store them in vec
%The first point of this segement start when reachs the length of (delta-res)
%The distance between the final point and end x2 is residue

% k -- the number of points that satisfy our condition
a=norm(x1-x2);

if (delta-res>a) || (a<1e-9)  
    residue=res+a;
    vec=[];
    k=0;
else
    k=1+floor((a+res-delta+1e-7)/delta);
    vec(:,1)=(delta-res)/a*(x2-x1)+x1;
    if k>1
    for i=1:k-1
    vec(:,i+1)=vec(:,i)+ (x2-x1)/a*delta;
    end
    
    residue=norm(x2-vec(:,k));
    else
     
    residue=norm(x2-vec(:,k));
    end

end
end

