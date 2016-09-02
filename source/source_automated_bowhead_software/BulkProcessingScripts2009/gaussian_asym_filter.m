%function h=gaussian_asym_filter(p2,p3,chc); 
    %p2: dimensions of NxM grid (vector)
    %p3: sigma values for x and y, respectively
    %chc:  if 1, vertical or horizontal output
    %    if greater than 2, diagonal output


function h=gaussian_asym_filter(p2,p3,chc); 


    siz   = (p2-1)/2;
    std   = p3;
     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     %arg   = -(x.*x + y.*y)/(2*std*std);
     
     if chc<=2
     arg   = -(x.*x)/(2*std(1).*std(1)) -(y.*y)/(2*std(2)*std(2));
     else
         arg   = -((x-y).^2)/(2*std(1).*std(1)) -((x+y).^2)/(2*std(2)*std(2));
    
     end


     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;