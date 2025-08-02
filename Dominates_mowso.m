function b=Dominates_mowso(x,y) %b=1 表示x数组中的值全都小于等于y数组中的值 且 至少有一个是小于y的 也就是说x支配y
    if isstruct(x)
        x=x.Cost;
    end
    if isstruct(y)
        y=y.Cost;
    end
     b=all(x<=y) && any(x<y);
end