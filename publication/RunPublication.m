function RunPublication
% Run all scripts needed to produce data and figure products for the publication describing
% P\lanetMag, DOI TBD.

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    GetNearFarMaxBdiff();
    CompareCoords();
    CalcBexc();
    CompareSatModels(0, 1, 16, 1);
    CompareSatModels(0, 1, 16, 0, 1);
    MakeHodograms();
    for planets=["Jupiter", "Saturn", "Uranus", "Neptune"]
        planet = char(planets);
        PrintSurfaceMap(planet, 0);
    end
end