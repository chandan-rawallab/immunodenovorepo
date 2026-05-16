import re

xml_path = "/home/amity/Documents/experiments/mqpar.xml"

with open(xml_path, 'r') as f:
    content = f.read()

# Lists of tags to truncate to 1 entry
tags_to_truncate = [
    "filePaths",
    "experiments",
    "fractions",
    "ptms",
    "paramGroupIndices",
    "referenceChannel"
]

for tag in tags_to_truncate:
    pattern = re.compile(rf'<{tag}>(.*?)</{tag}>', re.DOTALL)
    match = pattern.search(content)
    if match:
        inner_content = match.group(1)
        # Find all individual items (strings/shorts/booleans/ints)
        items = re.findall(r'<string>.*?</string>|<short>.*?</short>|<boolean>.*?</boolean>|<int>.*?</int>', inner_content)
        if items:
            # Keep only the first one
            new_inner = "\n      " + items[0] + "\n   "
            content = content.replace(inner_content, new_inner)

# Save result to a new temporary file to verify
with open("/home/amity/Documents/experiments/mqpar_single.xml", 'w') as f:
    f.write(content)
