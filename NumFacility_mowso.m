function pos = NumFacility_mowso(Q,VarMax,VarMin)
pos=rand(1,Q).*(VarMax-VarMin)+VarMin;
end