function [d,contrib,p] = BetaMixing(B,grp)
%Camden MacDowell - timeless
%uses one-way anovas to test whether betas are dominated by a single area,
%selectively mixed areas, or complete mix

%one-way anova
B=abs(B);
p = NaN(1,size(B,2));
contrib = NaN(1,size(B,2));
for i = 1:size(B,2)
   [p(i),~,stats] = anova1(B(:,i),grp,'off');
   c = multcompare(stats,'Display','off');
   idx = c(:,end)<0.05;
   if p(i)<0.05
       if sum(idx)==0
          contrib(i)=1;
       elseif numel(unique(c(idx,2)))==1 %only one is ever larger
           contrib(i)=3;
       else
           contrib(i)=2;
       end
   else
       contrib(i)=1;
   end     
end

%key | 1=none, 2=mixed, 3= solo
clear d; 
d(1) = sum(contrib==1)/numel(contrib);
d(2) = sum(contrib==2)/numel(contrib);
d(3) = sum(contrib==3)/numel(contrib);

end %function end


