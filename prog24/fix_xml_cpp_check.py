import sys
import xml.etree.ElementTree as ET

fname = sys.argv[1]
tree = ET.parse(fname)
root = tree.find("errors")
for node in root.findall('error'):
    if node.attrib["id"] == "missingIncludeSystem":
        root.remove(node)
tree.write(fname)