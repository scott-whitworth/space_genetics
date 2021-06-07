function series = angles(time,tFinal,coeff)
    series = coeff(1);
    for m = 1:(length(coeff)-1)/2
        series = series + coeff(2*m)*cos(2*pi*m*time/tFinal) + coeff(2*m+1)*sin(2*pi*m*time/tFinal);
    end
end

