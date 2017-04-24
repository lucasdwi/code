bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
for j = 1:size(bands,1)
    bandInd(j,1) = find(psdTrls.F>=bands{j,2}(1),1);
    bandInd(j,2) = find(psdTrls.F<=bands{j,2}(2),1,'last');
end
%%
for cond = 1:4
    for event = 1:2
        for anim = 1:sum(~cellfun(@isempty, rawPowArr(cond,:,event)))
            for chan = 1:4
                for trl = 1:size(rawPowArr{cond,anim,event},3)
                    thisTotPow = trapz(rawPowArr{cond,anim,event}(chan,bandInd(1,1):bandInd(end,end),trl));
%                     thisTotCoh = trapz(rawCohArr{cond,anim,event}(chan,bandInd(1,1):bandInd(end,end),trl));
                    for band = 1:6
                        thisBandPow = trapz(rawPowArr{cond,anim,event}(chan,bandInd(band,1):bandInd(band,2)));
%                         thisBandCoh = trapz(rawCohArr{cond,anim,event}(chan,bandInd(band,1):bandInd(band,2)));
                        relPow{cond,anim,event}(chan,band,trl) = thisBandPow/thisTotPow;
%                         relCoh{cond,anim,event}(chan,band,trl) = thisBandCoh/thisTotCoh;
                    end
                end
            end
        end
    end
end
for cond = 1:4
    for event = 1:2
        for anim = 1:sum(~cellfun(@isempty, rawCohArr(cond,:,event)))
            for chan = 1:6
                for trl = 1:size(rawCohArr{cond,anim,event},3)
%                     thisTotPow = trapz(rawPowArr{cond,anim,event}(chan,bandInd(1,1):bandInd(end,end),trl));
                    thisTotCoh = trapz(rawCohArr{cond,anim,event}(chan,bandInd(1,1):bandInd(end,end),trl));
                    for band = 1:6
%                         thisBandPow = trapz(rawPowArr{cond,anim,event}(chan,bandInd(band,1):bandInd(band,2)));
                        thisBandCoh = trapz(rawCohArr{cond,anim,event}(chan,bandInd(band,1):bandInd(band,2)));
%                         relPow{cond,anim,event}(chan,band,trl) = thisBandPow/thisTotPow;
                        relCoh{cond,anim,event}(chan,band,trl) = thisBandCoh/thisTotCoh;
                    end
                end
            end
        end
    end
end
%% Get average relPow and relCoh using same number of trials
nTrl = min(cellfun(@(x) size(x,3),relCoh),[],3);
for cond = 1:4
    for event = 1:2
        for anim = 1:sum(~cellfun(@isempty, rawCohArr(cond,:,event)))
            
        end
    end
end