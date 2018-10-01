def coding(base):
    if(base == 'A'):
        code = 0
    elif(base == 'C'):
        code = 1
    elif(base == 'G'):
        code = 2
    elif(base == 'T'):
        code = 3
    elif(base == '[C>A]'):
        code = 4
    elif(base == '[C>G]'):
        code = 5
    elif(base == '[C>T]'):
        code = 6
    elif(base == '[T>A]'):
        code = 7
    elif(base == '[T>C]'):
        code = 8
    elif(base == '[T>G]'):
        code = 9
    return code

def decoding(num):
    if(num == '0'):
        base = 'A'
    elif(num == '1'):
        base = 'C'
    elif(num == '2'):
        base = 'G'
    elif(num == '3'):
        base = 'T'
    elif(num == '4'):
        base = '[C>A]'
    elif(num == '5'):
        base = '[C>G]'
    elif(num == '6'):
        base = '[C>T]'
    elif(num == '7'):
        base = '[T>A]'
    elif(num == '8'):
        base = '[T>C]'
    elif(num == '9'):
        base = '[T>G]'
    return base

