function assignShapefunctions(obj)
for index=1:numel(obj)
    obj(index).shapeFunctions = lagrangeShapeFunctions(size(obj(index).edof,2), obj(index).numberOfGausspoints, obj(index).dimension);
end
end