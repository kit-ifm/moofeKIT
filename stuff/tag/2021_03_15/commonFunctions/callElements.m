function dataFE = callElements(obj,setupObject,flag)
dataFE = [];
for index1 = 1:obj.numberOfContinuumObjects
    continuumObj = evalin('base',obj.listContinuumObjects{index1});
     for index2 = 1:numel(continuumObj)
         if continuumObj(index2).flagNewton.(flag)
            elementName = strcat(continuumObj(index2).elementDisplacementType,continuumObj(index2).materialName,setupObject.integrator);
            dataFEindex = feval(elementName,continuumObj(index2),setupObject);
            dataFE = [dataFE, dataFEindex];
         end
     end
end