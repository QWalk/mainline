QW="../../bin/qwalk"
GOS="../../bin/gosling"

def check_errorbars(a,b,berr,sigma=3):
  if abs(a-b)/berr > sigma:
    return False
  return True

def check_sane(a,b,berr,sigma=3):
  if berr > 0.1:
    return False
  return True


