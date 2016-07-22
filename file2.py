import itertools


def patterntonumber(pattern):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=len(pattern))]
    return patterns.index(pattern)

print(patterntonumber('ATGCAA'))

def numbertopattern(number, k):
    patterns = [''.join(x) for x in itertools.product(['A', 'C', 'G', 'T'], repeat=k)]
    return patterns[number]

print(numbertopattern(5437, 7))
print(numbertopattern(5437, 8))

def computefrequency(text, k):
    result = [0 for i in range(4 ** k)]
    for i in range(len(text) - k + 1):
        result[patterntonumber(text[i:i+k])] += 1
    return result

print(computefrequency('ACGCGGCTCTGAAA', 2))
print(*computefrequency('aaa', 5), sep=' ')

def symbolToNumber(symbol):
    return {
        'A':0,
        'C':1,
        'G':2,
        'T':3
    }[symbol]

def patternToNumber2(pattern):
    if len(pattern) == 0:
        return 0
    prefix = pattern[:-1]
    symbol = pattern[-1:]
    return 4*patternToNumber2(prefix) + symbolToNumber(symbol)

print(patternToNumber2('TAATCGAAGATTAGCGCGA'))

def numberToPattern2(number,k):
    if k == 1:
        return numberToSymbol2(number)
    quotient, remainder = divmod(number,4)
    prefixPattern = numberToPattern2(quotient,k-1)
    symbol = numberToSymbol2(remainder)
    return prefixPattern+symbol

def numberToSymbol2(number):
    return {
        0:'A',
        1:'C',
        2:'G',
        3:'T'
    }[number]

print(numberToPattern2(7965,7))