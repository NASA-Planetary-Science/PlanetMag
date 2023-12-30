function RunPublication
    GetNearFarMaxBdiff();
    CompareCoords();
    CalcBexc();
    CompareSatModels(0, 1, 16);
    MakeHodograms();
end