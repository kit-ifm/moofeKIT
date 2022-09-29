def checkIntegerMatrix(input):
    if (input.isnumeric() and (sum(sum(input%1)))==0):
        print('input need to be of type integer even if it is classified as float')        
        
