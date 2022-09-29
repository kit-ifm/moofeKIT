function [patchData, meshData] = slimPatch(patchData,meshData)
%% remove double Vertices
[patchData.Vertices, indexm, indexn] = unique(round(patchData.Vertices*1e7)/1e7,'rows');
if isfield(patchData,'FaceVertexCData') && ~isempty(patchData.FaceVertexCData)
    patchData.FaceVertexCData = patchData.FaceVertexCData(indexm);
end
if size(patchData.Faces,1)==1
    patchData.Faces = indexn(patchData.Faces)';%does not work without transpose for single face (2D element)
else
    patchData.Faces = indexn(patchData.Faces);
end
if isfield(meshData,'Faces') && ~isempty(meshData.Faces)
    meshData.Faces = indexn(meshData.Faces);
end
%% remove double Faces
[patchData.Faces, indexm, indexn] = unique(patchData.Faces,'rows');
if isfield(meshData,'Faces') && ~isempty(meshData.Faces)
    [meshData.Faces, indexm, indexn] = unique(meshData.Faces,'rows');
end
end
