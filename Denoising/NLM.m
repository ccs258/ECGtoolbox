function [output]=NLM(input,t,f,h)
 % Input:
 %  input: signal to be filtered
 %  t: radio of search window
 %  f: radio of similarity window
 %  h: degree of filtering
 % Output:
 % output: filtered signal
 %  Author: Jose Vicente Manjon Herrera & Antoni Buades
 %  Date: 09-03-2006
 %
 %  Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
 %  "A non-local algorithm for image denoising"
 %
%%
 % length of the signal
 m = length(input);
 output =zeros(1,m);
 % Replicate the boundaries of the input image
  input1 = padarray(input,[1,t],'symmetric');
  input2 = input1(1,:);
 
 % Used kernel
 kernel = make_kernel(f);
 kernel = kernel / sum(sum(kernel));
 h=h*h;

 for j=1:m            

         j1 = j+ t;
         W1= input2(j1-f:j1+f);
         wmax=0; 
         average=0;
         sweight=0;
         B=t-f;
         smin = j1-B;
         smax = j1+B;   
         
         for s=smin:1:smax
                                               
%                   if(s==j) continue; end;
                                
                W2= input2(s-f:s+f);                
                 
                d = sum(sum(kernel.*(W1-W2).*(W1-W2)));
                w=exp(-d/h);                     
                if w>wmax                
                    wmax=w;                   
                end
                
                sweight = sweight + w;
                average = average + w*input2(s);                                  
         
         end
             
%         average = average + wmax*input2(j1);
%         sweight = sweight + wmax;
                   
        if sweight > 0
            output(1,j) = average / sweight;
        else
            output(1,j) = input(1,j);
        end                
 end

end

 
 
function [kernel] = make_kernel(f)   
kernel1=zeros(2*f+1,2*f+1); 
for d=1:f    
  value= 1 / (2*d+1)^2 ;    
  for i=-d:d
  for j=-d:d
    kernel1(f+1-i,f+1-j)= kernel1(f+1-i,f+1-j) + value ;
  end
  end
end
kernel1 = kernel1 ./ f;
kernel=sum(kernel1);
end