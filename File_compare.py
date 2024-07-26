
import filecmp
 
f1 = "Required_files/network1.txt"
f2 = "Required_files/network2.txt"
 
# deep comparison
result = filecmp.cmp(f1, f2, shallow=False)
print(result)