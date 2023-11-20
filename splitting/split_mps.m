function [mpData] = split_mps(mesh,mpData,directions)
to_split = zeros(length(mpData),1);
nmp = length(mpData);
nmp_new = nmp+1;
for mp=1:nmp 
    direction = directions(mp);
    if direction ~= 0
        mpData(nmp_new) = mpData(mp);
        position_disp = mpData(mp).lp(direction) * 0.5;
        new_lp = mpData(mp).lp * 0.5;
        new_lp0 = mpData(mp).lp0 * 0.5;
        new_vp = mpData(mp).vp * 0.5;
        new_vp0 = mpData(mp).vp0 * 0.5;
        new_mpM  = mpData(mp).mpM * 0.5;
        
        ids = [mp,nmp_new];
        normals = [-1,1];
        for i = 1:2
            cid = ids(i);
            normal = normals(i);
            mpData(cid).mpC(direction) = mpData(cid).mpC(direction) + (position_disp*normal);
            mpData(cid).lp(direction) = mpData(cid).lp(direction)*0.5;
            mpData(cid).lp0(direction)= mpData(cid).lp0(direction)*0.5;
            mpData(cid).vp = new_vp;
            mpData(cid).vp0 = new_vp0;
            mpData(cid).mpM = new_mpM;
        end
        nmp_new =  nmp_new + 1;
    end
end
