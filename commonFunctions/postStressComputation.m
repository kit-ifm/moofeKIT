function array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension)
plotObject = setupObject.plotObject;
stress = 0;
if strcmpi(plotObject.postPlotType, 'stress')
    stress = selectStress(stressTensor, setupObject, dimension);
elseif strcmpi(plotObject.postPlotType, 'D1') || strcmpi(plotObject.postPlotType, 'D2') || strcmpi(plotObject.postPlotType, 'D3')
    if strcmpi(plotObject.postPlotType, 'D1')
        indexD = 1;
    elseif strcmpi(plotObject.postPlotType, 'D2')
        indexD = 2;
    elseif strcmpi(plotObject.postPlotType, 'D3')
        indexD = 3;
    end
    stress = stressTensor.D(indexD);
end
array.Se = array.Se + N_k_I(k, :)' * stress * detJ * gaussWeight(k);
array.Me = array.Me + (N_k_I(k, :)' * N_k_I(k, :)) * detJ * gaussWeight(k);
end