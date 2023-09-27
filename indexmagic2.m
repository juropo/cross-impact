function [out] = indexmagic2(node,state,states)
%ITERATOR Summary of this function goes here
%   states should be a vector containing the number of states on each
%   variable
if(states(node)<state)
    error("Invalid state input. State exceeds the number of states for that node.")
end
if length(states)>1
    out = false(states);
else
    out = false(states,1);
end
call="out(";
for i=1:(node-1)
    call=call+":,";
end
call=call+state;
for i=(node+1):length(states)
    call=call+",:";
end
call=call+")=true;";
eval(call);
out=out(:);

end

