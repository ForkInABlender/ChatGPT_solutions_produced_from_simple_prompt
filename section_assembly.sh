# Dylan Kenneth Eliot & GPT-4-plugins (Beta Edition)

"""
Segment by section

* name
* then the data pertaining to that section

"""


BINARY="replace string with non-string path to binary"
echo $BINARY
sections=$(objdump -h $BINARY | grep '^[[:space:]]*[0-9]' | awk '{print $2}')

# Disassemble each section
for section in $sections; do
    echo "Disassembling section: $section"
    objdump -D -j $section $BINARY > "${section}.asm"
done
 
