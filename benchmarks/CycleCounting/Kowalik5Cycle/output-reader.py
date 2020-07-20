import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: output-reader.py [filename]")
    _, fn = sys.argv
    lines = []
    i = 2
    newline = []
    with open(fn, 'a') as f:
        while True:
            try:
                line = sys.stdin.readline()
                if line == None:
                    print("NONE")
                    break
                if line[:5] == "XXXXX":
                    if i < 2:
                        newline.append("TIMEOUT")
                    lines.append(",".join(newline) + "\n")
                    f.write(lines[-1])
                    newline = []
                    i = 0
                elif line[:5] == "#####":
                    i += 1
                    line = line.split(",")
                    newline.extend([x.split(":")[-1].strip().split("/")[-1] for x in line])
                    
            except:
                print("oops")
                break
