function [BxMoon, ByMoon, BzMoon] = RotateBspice(Bx, By, Bz, ets, spkParent, spkMoon)
    
    Bvec = [Bx; By; Bz];
    rotMat = cspice_pxform(spkParent, spkMoon, ets);
    BvecMat(:,1,:) = Bvec;
    BvecMoon = squeeze(pagemtimes(rotMat, BvecMat));
    BxMoon = BvecMoon(1,:);
    ByMoon = BvecMoon(2,:);
    BzMoon = BvecMoon(3,:);
end