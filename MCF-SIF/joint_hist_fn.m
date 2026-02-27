function [jhist_out]=joint_hist_fn(x1,x2)


[p,q]=size(x1);
jhist_out=zeros(256,256);
for ii=1:p
   for jj=1:q
      jhist_out(x1(ii,jj)+1,x2(ii,jj)+1)=jhist_out(x1(ii,jj)+1,x2(ii,jj)+1)+1;
   end
end