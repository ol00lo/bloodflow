import sys

fname=sys.argv[1]
print(fname)
fd=open(fname, "r")
lines=fd.readlines()
for i in range(len(lines)-2):
    line_cur = lines[i]
    line_next = lines[i+1]
    line_next_next=lines[i+2]
    line_cur_bare = line_cur.strip()
    line_next_bare = line_next.strip()
    if len(line_cur_bare)<4: 
        continue
    if line_cur_bare[:3] != "for":
        continue
    if line_cur_bare[-1] != ")":
        continue
    if line_next_bare[-1] == "{":
        continue
    if len(line_next_bare) < 4:
        continue
    num_of_for = 1
    for m in range(i+1, len(lines)-2):
        num_of_for += 1
        line_for = lines[m]
        lines_for_next = lines[m+1]
        lines_for_bare = line_for.strip()
        if(lines_for_bare[:3]== "for" and lines_for_next[-1] == "}"):
            break
    
    delta = len(line_next)-len(line_next_bare)-len(line_cur)+len(line_cur_bare)
    ifinal = i + 1
    sum_of_last = 0
    for k in range(i+1, len(lines)-2):
        ifinal += 1
        line_end = lines[k]
        line_end_next = lines[k + 1]
        line_end_bare = line_end.strip()
        line_end_next_bare = line_end_next.strip()
        if (len(line_end_bare.strip())!=0 and line_end_bare[-1]=="{"):
            sum_of_last+=1
        if (len(line_end_bare.strip())!=0 and line_end_bare[-1]=="}"):
            
            if sum_of_last > 1:
                sum_of_last-=1
            else:
                break

    for j in range(i+1, ifinal):
        lines[j] = lines[j][delta:]

fd=open(fname, "w")
for line in lines:
    fd.write(line)