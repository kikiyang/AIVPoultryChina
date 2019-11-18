function[G]=calculate_gravity(Mo,Md,distance,epsilon,beta,lambda)
[d11,d12]=size(Mo);
[~,d22]=size(Md);
if(d12~=1||d22~=1)  
    Mo=nanmean(Mo,2);
    Md=nanmean(Md,2);
end
k_exp = exp(- distance / (1000*lambda));
k_exp_norm = k_exp / (lambda^2);
if (epsilon==1&&beta==1)
    k_exp_norm=1 ./ ((distance/1000).^lambda);
end
for i=1:1:d11
     k_exp_norm(i,i)=0; %distance=0´¦£¬kernel(d)=0
end
index=1;
G=zeros(d11,d11);
while index<=d11
   G(index,:)=repmat(Mo(index).^epsilon,d11,1).*(Md.^beta).*k_exp_norm(:,index);
   index=index+1;
end