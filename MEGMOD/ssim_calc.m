function d = ssim_calc(x,y)

d=zeros(1,size(x,1));

for idx=1:size(x,1)
[~,lum_map]=ssim(x(idx,:)',y','Exponents',[1,0,0],'Radius',1.5,'DataFormat',"SS",'DynamicRange',1);
[~,con_map]=ssim(x(idx,:)',y','Exponents',[0,1,0],'Radius',1.5,'DataFormat',"SS",'DynamicRange',1);
[~,str_map]=ssim(x(idx,:)',y','Exponents',[0,0,1],'Radius',1.5,'DataFormat',"SS",'DynamicRange',1);
d(idx)=1-mean((lum_map).*(con_map).*(str_map));
end

end