




F = @(x) cos(x)/(1+exp(40*cos(x)));

asimpson(F,0,pi,10^(-7));






function x=asimpson(f,a,b,e)
s=0;
y=[a,b];
z=1;
k=0;

while length(y)>1
    for i=1:length(y)/2
        a1=y(2*i-1);
        b1=y(2*i);
        m=(a1+b1)/2;
        left=(m-a1)/6*(f(a1)+4*f((m+a1)/2)+f(m));
        right=(b1-m)/6*(f(m)+4*f((m+b1)/2)+f(b1));
        whole=(b1-a1)/6*(f(a1)+4*f(m)+f(b1));
        delta=left+right-whole;
        
        if abs(delta)<=15*e
            s=s+left+right+delta/15;
        else
            k=k+1;
            if k==1
                z=[a1,m,m,b1];
            else
                z=[z,a1,m,m,b1];
            end                
        end
    end
    y=z;    
    z=0;
    k=0;
end

x=s;
end