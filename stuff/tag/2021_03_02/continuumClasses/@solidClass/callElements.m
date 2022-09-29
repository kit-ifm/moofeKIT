function dataFE = callElements(obj,configObject)
dataFE = [];
for index = 1:numel(obj)
    elementName = strcat(lower(obj(index).elementDisplacementType),obj(index).materialName,configObject.integrator);
    dataFEindex = feval(elementName,obj(index),configObject);
    dataFE = [dataFE dataFEindex];
end
end
