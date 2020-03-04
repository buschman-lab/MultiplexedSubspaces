function [out_mat, norm_mat] = normalize_frobenius(in_mat,dim)
%normalize such that the frobenius norm 
assert(dim<=2,'dimesnion for normalization must be 1 or 2');

%column wise normalization
if dim ==1
    norm_mat=sqrt(sum(in_mat.^2,1));
    out_mat = bsxfun(@rdivide,in_mat,norm_mat);
else
    norm_mat=sqrt(sum(in_mat.^2,2));
    out_mat = bsxfun(@rdivide,in_mat,norm_mat);
end
        
end