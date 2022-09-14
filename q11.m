% for k=-100:1:100
%    w=(k/2)*(power(y,2)-power(x,2)) 
%    plot(w)
%    
%    
% end
[x,y] = meshgrid(-10:10);
u = y;
v = x;
l = streamslice(x,y,u,v);
axis tight