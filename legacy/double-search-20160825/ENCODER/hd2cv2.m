function [V,VI,C, row_weights, col_weights, H]=hd2cv2(HD,  M)
[b,c]=size(HD);
n=c*M;
r=b*M;
col_weights=zeros(1,n);
row_weights=zeros(1,r);
wcmax=max(sum(HD>=0));
wrmax=max(sum(HD'>=0));
C=zeros(n,wcmax);
V=zeros(r,wrmax);
VI=zeros(r,wrmax);
if nargout==6,
    H=zeros(r,n);
end; 

for i=1:b, 
    for j=1:c,
		if HD(i,j)>=0,
			for h=1:M	
                c_ind=(j-1)*M+mod((HD(i,j)+h-1),M)+1;  % column number
				r_ind=(i-1)*M+h;                 % row number      
				col_weights(c_ind)=col_weights(c_ind)+1;
				row_weights(r_ind)=row_weights(r_ind)+1;
                C(c_ind,col_weights(c_ind))=r_ind;
				V(r_ind,row_weights(r_ind))=c_ind;
				VI(r_ind,row_weights(r_ind))=col_weights(c_ind);
                if nargout==6, H(r_ind,c_ind)=1; end;
			end;
        end;
    end;
end;
	