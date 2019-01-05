function codeword=encode_QC_2BD(message, HB_U,HB_L, M)
% Generates random codeword for a LDPC code 
% with degree PCM in the form H= [P    U1 BD1 0
%                                 P1   P2     U2 BD2] where 
% P,P1,P2 are arbirary, U1,U2 are columns of weight 3 with 2 ones,
% BD1, BD2 are bidiagonal matrices

[bu,cu]=size(HB_U);
[bl,cl]=size(HB_L);
if cu~=cl;
    Z=zeros(bu,(cl-cu));
    Z=(Z==0).*(-1);
    HB_Uext=[HB_U Z];
end;

HB=[HB_Uext;HB_L]; %code parity-check matrix with upper part HB_U and lower part HB_L

if M==1,
    [~,~,~,~,~, H]=hd2cv2(HB,  M);
    g=g2h(H);
    codeword=mod(message*g,2);
    return;
end;

[b,c]=size(HB);
r=b*M;   % number of rows
n=c*M;   % number of columns
k=n-r;   % message length  

nu=cu*M;
ku=(cu-bu)*M; 

%%%%%Encode nu codesymbols corresponding to HB_U%%%%%%%%%%%%%
message_u=message(1:ku);
codeword(1:nu)=encode_QC(message_u,HB_U,M);
disp('first!');
%%%%%%%%%%%%%Encode n symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message_l=codeword(1:nu);
codeword(1:n)=encode_QC(message_l,HB_L,M);
disp('second');	
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



