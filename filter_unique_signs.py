

filename = "sign_file_diff.txt"
with open(filename) as sign_file:
    lines = sign_file.readlines()
    #unique lines
    lines = list(set(lines))

    pairs = [line.split() for line in lines]
    sums = [int(pair[0], 2) + int(pair[1], 2) for pair in pairs]
    to_be_removed = []
    for i in range(len(sums)-1):
        for j in range(i+1, len(sums)):
            if sums[i] == sums[j]:
                if (pairs[i][0] == pairs[j][1]) and (pairs[i][1] == pairs[j][0]):
                    to_be_removed.append(lines[j])

    unique_lines = [line for line in lines if line not in to_be_removed]
    with open("signs_unique_side.txt", "w") as unique_file:
        unique_file.writelines(unique_lines)
