
m=[1 2 3 4;5 6 7 8;9 10 11 12]
fcol=@(x)deal(x(:,1:2),x(:,3:4))
[a b]=fcol(m)

deal(1)
