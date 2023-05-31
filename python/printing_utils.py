def modify_printed_string(type,string):
    return "%s%s\033[0m"%(type,string)

def red(string):
    return modify_printed_string('\033[38;5;9m',string)

def green(string):
    return modify_printed_string('\033[38;5;2m',string)

def yellow(string):
    return modify_printed_string('\033[38;5;220m',string)

def orange(string):
    return modify_printed_string('\033[38;5;208m',string)

def blue(string):
    return modify_printed_string('\033[38;5;12m',string)

def magenta(string):
    return modify_printed_string('\033[38;5;5m',string)

def cyan(string):
    return modify_printed_string('\033[38;5;6m',string)

def bold(string):
    return modify_printed_string('\033[1m',string)

def prettydict(d, indent=4, color=blue):
    space = max([0]+[len(str(x)) for x in d])+2
    print('')
    for key, value in d.items():
        print(color(" "*indent + str(key)),end='')
        if isinstance(value, dict):
            prettydict(value, indent=len(" "*(space)), color=color)
        else:
            print(color(" "*(space-len(key)) + str(value)))
