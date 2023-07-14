function s = jmnum2str(n,digit)
s = num2str(n);
for i=1:digit
    j=digit-i;
    if n<10^(j)
        s = ['0' s];
    end
end
        