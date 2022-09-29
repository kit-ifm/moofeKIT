function dataFE = callElements(obj,setupObject,flag)
dataFE = [];
for index1 = 1:obj.numberOfContinuumObjects
    continuumObject = obj.listContinuumObjects{index1};
    if continuumObject.flagNewton.(flag)
        elementName = strcat(continuumObject.elementDisplacementType,continuumObject.materialName,setupObject.integrator);
        dataFEindex = feval(elementName,continuumObject,setupObject);
        dataFE = [dataFE, dataFEindex];
    end
end