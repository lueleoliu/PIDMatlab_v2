function score = grading(ts,tr,mp,osc)
%This function is defined to return the score of a group of parameters in Pareto matrix
score = 0;

if ts<15
    score = score + 5;
elseif ts<20
    score = score + 4;
elseif ts<25
    score = score + 3;
elseif ts<30
    score = score + 2;
elseif ts<35
    score = score + 1;
else
    score = score + 0;
end
    
if tr<5
    score = score + 5;
elseif tr<6
    score = score + 4;
elseif tr<7
    score = score + 3;
elseif tr<8
    score = score + 2;
elseif tr<9
    score = score + 1;
else
    score = score + 0;
end

if mp<10
    score = score + 0;
elseif mp<20
    score = score + 5;
elseif mp<25
    score = score + 4;
elseif mp<30
    score = score + 3;
elseif mp<35
    score = score + 2;
elseif mp<40
    score = score + 1;
else
    score = score + 0;
end

if osc < 1
    score = score + 0;
elseif osc<3
    score = score + 5;
elseif osc<4
    score = score + 4;
elseif osc<5
    score = score + 3;
elseif osc<6
    score = score + 2;
elseif osc<7
    score = score + 1;
else 
    score = score + 0;
end

end


