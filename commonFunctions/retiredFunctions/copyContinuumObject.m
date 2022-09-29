function [continuumObjectOut] = copyContinuumObject(dofObject,continuumObjectIn)
if ~isempty(dofObject)
    a = [class(continuumObjectIn), '(dofObject,', '''nolisting''',')'];
else
    a = class(continuumObjectIn);
end
continuumObjectOut = eval(a);  %create default object of the same class as a. one valid use of eval
for p = properties(continuumObjectIn).'  %copy all public properties
    propertyName = p{1};
    if isobject(continuumObjectIn.(propertyName))
        continuumObjectOut.(propertyName) = copyContinuumObject([],continuumObjectIn.(propertyName));
    else
        try %may fail if property is read-only
            continuumObjectOut.(propertyName) = continuumObjectIn.(propertyName);
        catch
            warning('failed to copy property: %s', propertyName);
        end
    end
end
end