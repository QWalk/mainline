from __future__ import print_function
import csv
import math

QW="../../bin/qwalk"
GOS="../../bin/gosling"

def check_errorbars(a,aerr,b,berr,sigma=3):
  if abs(a-b)/math.sqrt(aerr**2+berr**2) > sigma:
    return False
  return True

def check_sane(berr):
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

def compare_result_ref(ref_data,dat_properties,sigmas):
  success={}
  for k in ref_data.keys():
    success[k]=check_errorbars(ref_data[k][0],
                               ref_data[k][1],
                               dat_properties[k]['value'][0],
                               dat_properties[k]['error'][0],
                               sigmas[k])
    success[k+'sane']=check_sane(dat_properties[k]['error'][0])
  return success

def summarize_results(ref_data,dat_properties,success,systems,methods,descriptions):
  reports=[]
  for k in ref_data.keys():
    report={}
    report['system']=systems[k]
    report['method']=methods[k]
    report['description']=descriptions[k]
    report['quantity']=k
    report['result']=dat_properties[k]['value'][0]
    report['error']=dat_properties[k]['error'][0]
    report['reference']=ref_data[k][0]
    report['err_ref']=ref_data[k][1]
    report['passed']=success[k] and success[k+'sane']
    reports.append(report)
  return reports
