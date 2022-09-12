function [cs, sn] = GivensRotationCoeff(v1, v2)
  if (v1==0)
    cs = 0;
    sn = 1;
  else
    t = abs(v1) * sqrt(1 + (v2/v1)^2);
    cs = abs(v1) / t;
    sn = cs * v2 / v1;
  end
end