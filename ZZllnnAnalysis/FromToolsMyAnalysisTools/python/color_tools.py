
# ------------------------------------------------------------------------------------
# colors for printout

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[33m'
    WARNING_H = '\033[33m'
    WARNING_L = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    @staticmethod
    def disable():
        bcolors.HEADER = ''
        bcolors.OKBLUE = ''
        bcolors.OKGREEN = ''
        bcolors.WARNING = ''
        bcolors.WARNING_H = ''
        bcolors.WARNING_L = ''
        bcolors.FAIL = ''
        bcolors.ENDC = ''

    @staticmethod
    def enable():
        bcolors.HEADER = '\033[95m'
        bcolors.OKBLUE = '\033[94m'
        bcolors.OKGREEN = '\033[92m'
        bcolors.WARNING = '\033[33m'
        bcolors.WARNING_H = '\033[33m'
        bcolors.WARNING_L = '\033[93m'
        bcolors.FAIL = '\033[91m'
        bcolors.ENDC = '\033[0m'

def blue(string):
    return bcolors.OKBLUE + string + bcolors.ENDC

def ok(string):
    return bcolors.OKGREEN + string + bcolors.ENDC

def warning(string):
        return bcolors.WARNING + string + bcolors.ENDC

def warningH(string):
        return bcolors.WARNING_H + string + bcolors.ENDC

def warningL(string):
        return bcolors.WARNING_L + string + bcolors.ENDC



def error(string):
        return bcolors.FAIL + string + bcolors.ENDC
