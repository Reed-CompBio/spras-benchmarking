
import filecmp

f1 = r"C:\Users\dnachreiner\Desktop\BDS_Project\Scripts\Raw_data\network1.txt"
f2 = r"C:\Users\dnachreiner\Desktop\BDS_Project\Scripts\Raw_data\network2.txt"

# deep comparison
result = filecmp.cmp(f1, f2, shallow=False)
print(result)