function txtOUT = jjnum2str(valIN,n)
%Convert values to txt up to .00 decimal point
%Jangjin Kim, 2011-May-31

txtOUT = num2str(round(valIN * 10^n) / 10^n);