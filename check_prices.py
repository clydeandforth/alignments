#!/scratch/b.bss81c/anaconda3/bin/python
# which price is better value for money tri_shipping, drone or premimum

import sys

value = sys.argv[1]

def tri_shipping(value):
  num = 20
  if int(value) <= 2:
    return float(value) * 1.5 + num
  elif int(value) > 2 and int(value) <= 6:
    return float(value) * 3.00 + num
  elif int(value) > 6 and int(value) <=10 :
    return int(value) * 4 + num
  else:
    return float(value) * 4.75 + num
  
def drone(value):
   if int(value) <= 2:
    return float(value) * 4.5 
   elif int(value) > 2 and int(value) <= 6:
    return float(value) * 9.00 
   elif int(value) > 6 and int(value) <=10 :
    return int(value) * 12 
   else:
    return float(value) * 14.75 

  

def overall():
  value1=tri_shipping(value)
  value2=drone(value)
  premium=125
  if value1 > premium and value2 > premium:
    print("use premium which costs $" + str(premium))
  elif value1 > value2:
    print("use drone which costs $" + str(value2))
  elif value1 < value2:
    print("use ground which costs $"+str(value1)) 
  else:
    print("use any")
    
def main():
  overall()
main()    
  

