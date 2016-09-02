%%Illustrates how variables are passed in nested functions.
%% Key point, variables defined in nested functions can be passed to higher
%% levels, if upper level refers to variable in question.  Otherwise, not
%% passed to workspace.  Kinda looser.

function varScope1
x = 5;
y=nestfun1;

disp(y)
disp(z)
disp(zz);
%disp(xx);


keyboard
   function y=nestfun1
      y=nestfun2;
      xx=zeros(10,1);
        zz=10;
        
      function y=nestfun2
         y = x + 1;
         z=4;
        
      end
   end


end