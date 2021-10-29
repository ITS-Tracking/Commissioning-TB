import re

with open('data_all.lst') as f:
    lines = f.read().splitlines()
    for i in range(0,len(lines)):
        fname = lines[i][0:-5]
        num = 0 
        if fname[-2] == '_':
            num = int(fname[-1])
            fname = lines[i][0:-7]    
        r = re.compile(fname+"\_[1-9].root")
        others = list(filter(r.match, lines[i+1:]))
        if len(others) == 0:
            print(lines[i])
        else:
            printalo = True
            for other in others:
                otherNum = int(other[-6])
                if otherNum > num:
                    printalo = False
            if printalo:
                print(lines[i])

