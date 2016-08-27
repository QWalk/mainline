from __future__ import print_function
import csv

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

def print_results(reports):
  for report in reports:
    print("%10s %15s %5s %6r %f+/-%f %f+/-%f"%(report['method'],report['quantity'],report['system'],report['passed'],report['result'],report['error'],report['reference'],report['err_ref']))

def save_results(reports):
  fieldnames = ['method','quantity','system','description','passed','result','error','reference','err_ref']
  with open('report.csv','w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writerow(dict(zip(fieldnames,fieldnames)))
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writerows(reports)

