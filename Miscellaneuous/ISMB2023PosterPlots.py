
# take the initial items from the paperFigure01.py

for accession in accessionList:
    annotation = allMeta[accession]
    if annotation['tissueInfo'][0] == 'pancreas':
        print(accession)


# print tissue names with order

tissueNameList = []
for accession in accessionList:
    annotation = allMeta[accession]
    tissueNameList.append(annotation['tissueInfo'][0])

book = sorted((name, index) for index,name in enumerate(tissueNameList))
for i, name in enumerate(book):
    print(i, name, accessionList[name[1]])

# in this sorted list:
'''
2 of 1
1 of 3
1 of 1
6 of 10 (close family)
2 of 12 (close family)
13-31: 1 of each (cell line)
2 of 33: (close family)
2 of 35: close family
3 of 39
40-46: 1 of each (cell line)
3 of 49
50-51: one of each
3 of 54: close family
1 of 55
2 of 57
2 of 59
1 of 60 (cardiac muscle)
2 of 62
4 of 66
1 of 67
3 of 70
2 of 72
1 of 73 (cortex brain)
74-77: one of each
1 of 78
2 of 80
8 of 89 (tissue, not cell type)
3 of 92
5 of 97
2 of 101
1 of 102
20 of 122 (heart ventricle)
1 of 123
4 of 127 (close)
128-130: 1 of each
3 of 133
2 of 135, 4 of 140-143 (lung)
1 of 136
3 of 139 (liver)
144-145: 1 of each
2 of 147
4 of 151 (muscle)
1 of 152
3 of 155 (family of 6-10)
1 of 156 (family of above)
3 of 158 (neuron)
159-162: one of each
5 of 167
1 of 168
2 of 170
4 of 174 (close family)
1 of 175
2 of 177 (muscle)
2 of 179
3 of 182 (liver)
3 of 187
2 of 189 (muscle)
2 of 191
1 of 192
1 of 193
8 of 201 (spleen)
5 of 206 (stomach)
207-209: one of each
2 of 211 (brain)
212-213: one of each
2 of 215
3 of 218
1 of 219
4 of 223
224-225: one of each
5 of 230
2 of 232
1 of 233
'''



