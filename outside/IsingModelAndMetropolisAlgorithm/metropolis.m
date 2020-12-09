function spin = metropolis(spin, kT, J)
%METROPOLIS The Metropolis algorithm.
%   spin = METROPOLIS(spin, kT, J) runs the Metropolis algorithm on a
%   configuration of spins with an coupling coefficient |J| at a
%   temperature |kT|. |spin| is a matrix of +/- 1's.

%   Copyright 2017 The MathWorks, Inc.

numIters = 2^8 * numel(spin);
h = figure;
pcolor(spin)
colormap(gray)
filename = 'test.gif';
c = 1;
for iter = 1 : numIters
    
    % Pick a random spin
    linearIndex = randi(numel(spin));
    [row, col]  = ind2sub(size(spin), linearIndex);
    
    % Find its nearest neighbors
    above = mod(row - 1 - 1, size(spin,1)) + 1;
    below = mod(row + 1 - 1, size(spin,1)) + 1;
    left  = mod(col - 1 - 1, size(spin,2)) + 1;
    right = mod(col + 1 - 1, size(spin,2)) + 1;
    
    neighbors = [      spin(above,col);
        spin(row,left);                spin(row,right);
                       spin(below,col)];
    
    % Calculate energy change if this spin is flipped
    dE = 2 * J * spin(row, col) * sum(neighbors);
    
    % Boltzmann probability of flipping
    prob = exp(-dE / kT);
    
    % Spin flip condition
    if dE <= 0 || rand() <= prob
        spin(row, col) = - spin(row, col);
    end
    c = c+1;
    if iter == 1 | c == 100
        pcolor(spin);
        colormap(gray)
        drawnow
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if iter == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','writemode','append');
        end
        c = 1;
    end
    
end
% figure
% pcolor(spin)
% colormap(gray)
