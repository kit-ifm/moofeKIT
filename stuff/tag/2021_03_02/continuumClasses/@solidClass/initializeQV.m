function initializeQV(obj)
for index = 1:numel(obj)
    if isempty(obj(index).QN)
        obj(index).QN = obj(index).QR;
    end
    if isempty(obj(index).QN1)
        obj(index).QN1 = obj(index).QR;
    end
    if isempty(obj(index).VN)
        obj(index).VN = zeros(size(obj(index).nodes,1), obj(index).dimension);
    end
end
end