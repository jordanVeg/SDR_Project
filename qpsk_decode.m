function dec = qpsk_decode(bit0, bit1)
    if bit0 == -1 && bit1 == -1
        dec = -3;
    elseif bit0 == -1 && bit1 == 1
        dec = -1;
    elseif bit0 == 1 && bit1 == -1
        dec =1;
    else
        dec=3;
    end
end