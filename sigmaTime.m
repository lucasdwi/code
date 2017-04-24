function [sigmaT] = sigmaTime(gam,bet)
sigSq = normConst(gam,bet)^2/nZero(gam,bet)*(bet^2*(nZero(gam,bet-1)/normConst(gam,bet-1)^2)+gam^2*(nZero(gam,bet-1+gam)/normConst(gam,bet-1+gam)^2)-2*gam*bet*(nZero(gam,bet-1+gam/2)/normConst(gam,bet-1+gam/2)^2));
sigmaT = sqrt(sigSq);
    function [alp] = normConst(gam,bet)
        alp = 2*((exp(1)*gam)/bet)^(bet/gam);
    end
    function [M] = mZero(gam,bet)
        M = (normConst(gam,bet)/(2*pi*gam))*gamma((bet+1)/gam);
    end
    function [N] = nZero(gam,bet)
        N = (2/(2^(1/gam)))*mZero(gam,2*bet);
    end
%     function [x] = gammaHat(z)
%         x = gamma(z)/2^z;
%     end
end