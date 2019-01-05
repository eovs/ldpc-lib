function codeword=encode_QC(message, HB, M)
% Generates random codeword for a LDPC code 
% with degree PCM in the form H= [P | U | BD] where 
% P is arbirary, U is column of weight 3 with 2 ones,
% BD is bidiagonal 
if M==1,
    [~,~,~,~,~, H]=hd2cv2(HB,  M);
    g=g2h(H);
    codeword=mod(message*g,2);
    return;
end;
[b,c]=size(HB);
r=b*M ;  % number of rows
n=c*M ;  % number of columns
k=n-r ;  % message length 

% find positive degrees in (c-b+1)st column
L=zeros(3,2);
codeword=zeros(1,n);
j=0;
for i=1:b 
    if HB(i,c-b+1)>=0,
        j=j+1;
        if j>3, error('weight of (c-b+1)st column is too large'); end;
        L(j,1)=i;
        L(j,2)=HB(i,c-b+1);
    end;
end;
if j<3, error('weight of (c-b+1)st column is too small'); end;
if L(1,2)~=0 || L(3,2)~=0 || L(2,2)==0,  
    error('(c-b+1)st column must be [0...A...0]');
end;
% Generate info part of codeword
%message=randi([0 1], 1,k);
%load message;

codeword(1:k)=message; 
% Compute partial syndrome
synd=zeros(1,r);
for j=1:c-b    % for all columns 
    y=codeword((j-1)*M+(1:M)); % read block
    for i=1:b  % for all rows
        % circulate and add
        if HB(i,j)>=0
            yc=cyclic_shift_left(y,HB(i,j));
            synd((i-1)*M+(1:M))=synd((i-1)*M+(1:M))+yc;
        end;
    end;
end;
synd=mod(synd,2);
% Compute sum of syndrom components
sumsynd=mod(sum(reshape(synd,M,b),2),2)';
% One check block is known
codeword(k+(1:M))=cyclic_shift_left(sumsynd,M-HB(L(2,1),c-b+1));
% Partial syndrom modification
%synd((L(1,1)-1)*M+(1:M))=synd((L(1,1)-1)*M+(1:M))+codeword(k+(1:M));
synd((L(2,1)-1)*M+(1:M))=mod(synd((L(2,1)-1)*M+(1:M))+sumsynd,2);
synd((L(3,1)-1)*M+(1:M))=mod(synd((L(3,1)-1)*M+(1:M))+codeword(k+(1:M)),2);
% recursion
for i=c-b+2:c
    j=i-c+b-2;
	codeword((i-1)*M+(1:M))=mod(synd(j*M+(1:M))+codeword((i-2)*M+(1:M)),2);
end;	
% CHECK for being a proper codeword
synd=zeros(1,r);
for j=1:c,
    y=codeword((j-1)*M+(1:M)); % read block
    for i=1:b,
        if HB(i,j)>=0
            yc=cyclic_shift_left(y,HB(i,j));
            synd((i-1)*M+(1:M))=synd((i-1)*M+(1:M))+yc;
        end;
    end; 
end;
synd=mod(synd,2);
if ~all(synd==0), error('bad coding'); end;

% --------------------------------------------------%
function y=cyclic_shift_left(x,s)
% y is a cyclic shift of x by s positions left
M=length(x);
y=[x(s+1:M) x(1:s)];


