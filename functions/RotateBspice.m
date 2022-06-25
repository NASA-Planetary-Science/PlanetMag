function [BxMoon, ByMoon, BzMoon] = RotateBspice(Bx, By, Bz, ets, spkParent, spkMoon)
    
    Bvec = [Bx; By; Bz];
    npts = length(ets);
    rotMat = cspice_pxform(spkParent, spkMoon, ets);
    BvecMoon = zeros(3,npts);
    parfor i=1:npts
        BvecMoon(:,i) = rotMat(:,:,i) * Bvec(:,i);
    end
    BxMoon = BvecMoon(1,:);
    ByMoon = BvecMoon(2,:);
    BzMoon = BvecMoon(3,:);
end