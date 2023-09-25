%Ellipse of Insignificance MATLAB function. DR Grimes 2023
%Based on work in: The ellipse of insignificance, a refined fragility index for ascertaining 
%robustness of results in dichotomous outcome trials

%https://doi.org/10.7554/eLife.79573

%Usage: af = events experimental arm, bf = events control arm
%cf = total experimental arm, df = total control arm;
%alpha = significance level


function y = checksig(af,bf,cf,df,alpha)
   a = af;
   b = cf - af;
   c = bf;
   d = df - bf;
   
   
   chival = @(v) chi2cdf(v,1) - (1-alpha);
v = fzero(chival,0);
n = a + b + +c +d; 
tstat_test = (n.*((a.*d - b.*c).^2))./((a+b).*(c+d).*(a+c).*(b+d)) - v;
pvalchk = chi2cdf(tstat_test + v,1);
hr = (a/(a+b))/(c/(c+d));
lower = exp(log(hr) - 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)));
uppr = exp(log(hr) + 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)));


 if tstat_test >= 0
        [xi yi dmin qvec1 qvec2 qvec3] = eoi(af,bf,cf,df);
    else
        'Not significant'; % Set x to an empty array or handle the case as needed.
        assignin('base', 'val', [1-pvalchk 0 0 0 0 0 0 hr lower uppr]);
    end


end

function [xi yi dmin qvec1 qvec2 qvec3] = eoi(af,bf,cf,df)
   a = af;
   b = cf - af;
   c = bf;
   d = df - bf;
   alpha = 0.05;  

n = a + b + +c +d;
   chival = @(v) chi2cdf(v,1) - (1-alpha);
v = fzero(chival,0);
tstat_test = (n.*((a.*d - b.*c).^2))./((a+b).*(c+d).*(a+c).*(b+d)) - v;
pvalchk = chi2cdf(tstat_test + v,1);


       %ellipse terms
Ae = (c+d).*((c+d).*n + (a+b).*v);
Be = 2*(a+b)*(c+d)*(n-v); 
Ce = (a+b).*((a+b).*n + (c+d).*v);
De = (c+d).*(2.*(b*c-a*d).*n + (a+b).*(b-a + d -c).*v);
Ee = (a+b).*(2.*(b*c-a*d).*n + (c+d).*(a-b + c -d).*v);
Fe = ((b*c-a*d)^2).*n - (a+b)*(a+c)*(b+d)*(c+d)*v; 

%find FECKUP vector points
syms xp yp
eq1=(2*Ae.*xp + Be.*yp +De).*yp - xp.*(Be.*xp + 2*Ce.*yp + Ee);
eq2=Ae*(xp.^2) + Be.*xp.*yp + Ce.*(yp.^2) + De.*xp + Ee.*yp + Fe; 
eqs = [eq1, eq2];
[xp,yp]=vpasolve(eqs,[xp,yp]);
xp = double(xp);
yp = double(yp);
dv = sqrt((xp.^2) + (yp.^2)); 
i = find(dv == min(dv)); 
xp = xp(i);
yp = yp(i); 
FCKP = dv(i);
clear eq1 eq2 eqs dv i 

%find xi point - simple quadratic form
A1 = (c+d).*((c+d).*n + (a+b).*v);
B1 = (c+d).*(2.*(b*c-a*d).*n + (a+b).*(b-a + d -c).*v);
C1 = ((b*c-a*d)^2).*n - (a+b)*(a+c)*(b+d)*(c+d)*v; 
X1 = (-B1 + sqrt(B1.^2 - 4*A1*C1))./(2*A1); 
X2 = (-B1 - sqrt(B1.^2 - 4*A1*C1))./(2*A1);
g = [X1 X2];
gc = abs(g);
xi = g(find(gc == min(gc))); 
clear A1 B1 C1 X1 X2 g

%find yi point - simple quadratic form
A2 = (a+b)*((a+b)*n + (c+d)*v);
B2 = (a+b)*(2*(b*c-a*d)*n + (a-b+c-d)*(c+d)*v);
C2 = ((b*c-a*d)^2).*n - (a+b)*(a+c)*(b+d)*(c+d)*v; 
Y1 = (-B2 + sqrt(B2.^2 - 4*A2*C2))./(2*A2); 
Y2 = (-B2 - sqrt(B2.^2 - 4*A2*C2))./(2*A2);
g2 = [Y1 Y2];
g2c = abs(g2);
yi = g2(find(g2c == min(g2c))); 
clear A2 B2 C2 Y1 Y2 g2

%Find all relevant terms

dmin = floor(abs(xp) + abs(yp));
qvec = [(1 - (a+b - abs(xi))./(a+b)) (1 - (c+d - abs(yi))./(c+d)) (1 - (n - abs(dmin))./(n))];
qvec1 = qvec(1);
qvec2 = qvec(2);
qvec3 = qvec(3); 

hr = (a/(a+b))/(c/(c+d)); %is rr 

%https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html#:~:text=Therefore%2C%20computing%20the%20confidence%20interval,confidence%20interval%20for%20the%20RR.

lower = exp(log(hr) - 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)));
uppr = exp(log(hr) + 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)));


vals = [1-pvalchk xi yi dmin 100*qvec1 100*qvec2 100*qvec3 hr lower uppr];

assignin('base', 'val', vals);
       
   
end

