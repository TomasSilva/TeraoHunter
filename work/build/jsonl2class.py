import json, sys
fn, outprefix = sys.argv[1], sys.argv[2]
for i, l in enumerate(open(fn)):
    r = json.loads(l)
    with open(f"{outprefix}{i}.txt", "w") as f:
        f.write("sizes " + " ".join(map(str, r["sizes"])) + "\n")
        for B in r["blocks"]: f.write("point " + " ".join(map(str, B)) + "\n")
