function fraction=laser_reduction(x)

para = [ 7.63395370e-02, -6.58170082e-01,  2.13538100e+00, -3.13548207e+00,  1.95182504e+00, -5.38481103e-01,  3.67228491e-01, -3.07200121e-04];
f=@(x) polyval(para,x*1.1)/0.2084; % The particle radius was 1.1, and the max of the curve was 0.2084
if x<2.1
    fraction=f(x);
else 
    fraction=0;
end
end