function codeword=encode_lowtr(message, HB, M)
% Generates random codeword for a LDPC code 
% with degree PCM in the form H= [P | LT] where 
% P is arbirary, LT is low-triangular 
if M==1,
    [~,~,~,~,~, H]=hd2cv2(HB,  M);
    g=g2h(H);
    codeword=mod(message*g,2);
    return;
end;
[b,c]=size(HB);
r=b*M   % number of rows
n=c*M   % number of columns
k=n-r   % message length  

codeword=zeros(1,n);

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
synd
pause
synd=mod(synd,2);

% One check block is known
codeword(k+(1:M))=synd(1:M);
tmp=synd(1:M);
%recursion
j=2;
for i=c-b+2:c
    %j=i-c+b;
    %%%Partial syndrome modification%%%%%
    %%%Compute a new syndrome for the current block of the codeword%%%%%%
    for t=j:b  % for all rows
        % circulate and add
        if HB(t,i-1)>=0
            yc=cyclic_shift_left(tmp,HB(t,i-1));
            synd((t-1)*M+(1:M))=synd((t-1)*M+(1:M))+yc;
        end;
    end;
    synd=mod(synd,2)
    
    tmp=synd((j-1)*M+(1:M))
    pause
	codeword(k+(j-1)*M+(1:M))=tmp;
    j=j+1;
end;
%return;
  if 1
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
end;
% --------------------------------------------------%
function y=cyclic_shift_left(x,s)
% y is a cyclic shift of x by s positions left
M=length(x);
y=[x(s+1:M) x(1:s)];


