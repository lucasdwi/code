permlist = perms(1:3);
for ei = 1:2
   for ti = 1:length(trls{1,ei}.trial)
      [pe(ei,ti),hist(ti,:,ei),z] = pec(trls{1,ei}.trial{1,ti}(1,:),3,1);
      for ii = 1:size(z,1)
          for jj=1:length(permlist)
              if (abs(permlist(jj,:)-z(ii,:)))==0
                  pat(ti,ii,ei) = jj;
              end
          end
      end
   end
end