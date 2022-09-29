function checkIntegerMatrix(input)
assert(isnumeric(input) & (sum(sum(mod(input,1))) == 0),'input need to be of type integer even if it is classified as float');
end